#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>

#include "ksw.h"
#include "bwa.h"
#include "bwt.h"
#include "aln.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"


/* calculate  maximun gap length */

static inline int cal_max_gap (const opt_t *opt , int len)
{	
	int  l_del = (int)((double)(len*opt->a - opt->o_del)/opt->e_del + 1.);
	int  l_ins = (int)((double)(len*opt->a - opt->o_ins)/opt->e_ins + 1.);
	int  l =  l_del  >  l_ins  ?  l_del : l_ins ;
	l = l > 1 ?  l : 1 ;
	return  l < opt->w << 1 ? l : opt->w << 1;
}

/*  
 *  cm_meraln : (len -  k  + 1) k-mer seeds align  to generate  result.
 *
 */

int   cm_meraln( aln_seed_v *kv_seed ,  const seq_t *pseq , const opt_t *opt)
{
	unsigned char *bseq ;

	int	i  =  0  ;
	
	aln_seed_t  s_aln;
	
	char2nt4(&bseq,pseq->len,pseq->seq,0);
	seq_complement(pseq->len ,bseq);
	
	for( i  = 0 ;  i  <  pseq->len - opt->l_seed + 1 ; i++){
		bwtint_t  k , l;
		k = 0 , l = opt->fr->idx->bwt->seq_len ;
		/* 
		 * 	bwt_match_exact_cm function similar to bwt_match_exact in bwa.c ,but   sequence (bseq) is  complenment.
		 * 	generate k-mer seeds  to call the bwt_match_exact_cm
		 */
		if(bwt_match_exact_cm( opt->fr->idx->bwt , i , opt->l_seed ,bseq  , &k , &l)){ // set  MAX_INTERVAL  , to filter redundancy k-mer hit..
			bwtint_t  j ,pos;
			for( j = k ;  j < l + 1 ; j++){
				/*
				 *  This section very slow ,  acquaire some coordinate by method of iteration(SA =? iSA).
				 */
				pos = bwt_sa(opt->fr->idx->bwt,j);
				s_aln.ref_b =  pos , s_aln.ref_e = pos + opt->l_seed ;
				s_aln.query_b =  pseq->len - i - opt->l_seed , s_aln.query_e =   pseq->len - i ;
				kv_push(aln_seed_t,*kv_seed,s_aln);
			}
		}
	}
	free(bseq);
	return  kv_size(*kv_seed);
}


/*
 *   this is bucket sort  from  ksort.h
*/

#define flt_fuc(a,b) ((a).ref_b < (b).ref_b)
KSORT_INIT(cm_seed_flt,aln_seed_t,flt_fuc)

/*
 *     @abstract: merge  adjacent coordinate seed to chain
 */

void   cm_mergechain(aln_chain_v *av ,  aln_seed_v  *kv_seed)
{
	int  j = 0  , i; 
	
	aln_seed_t  *s_aln;
	bwtint_t    max =  0 ;
//       first , regular merge seed .
	ks_introsort(cm_seed_flt,kv_size(*kv_seed),kv_seed->a);
	s_aln = kv_seed->a ;
	max =  s_aln->ref_e ; 
	push_chain_seed(*av,*s_aln,aln_seed_t,0);

	for ( i = 1 ;  i < kv_size(*kv_seed) ; i++){
		aln_chain_t  *p  =  av->a + av->n ; 
		s_aln = kv_seed->a + i ;
		
		for(  j = 0   ;  j  < p->n ; j++){
	/*
	*	  check whether chain have unmerge seed.
	*	  eg:  ref: ------ACCCTGATCCCTGA------
	*	  	reads     ACCCTGATCCCTGACCCTGATCCC
	*	 example:   tandemCNV 
	*/ 
			aln_seed_t *p_aln =  p->a + j;
			if( max < p_aln->ref_e )  max = p_aln->ref_e ;
			if( p_aln->ref_e <= s_aln->ref_b ) continue ;

			if( s_aln->query_b - p_aln->query_b == s_aln->ref_b - p_aln->ref_b){
				p_aln->ref_e =  s_aln->ref_e ;
				p_aln->query_e = s_aln->query_e ;
				break;
			}
		}
		if( j == p->n ){
			//    Two adjacent  seeds  offset  < av->offset   , those belong to a chain .  
			if( max + av->offset > s_aln->ref_b ){
				s_aln->len =  s_aln->query_e - s_aln->query_b ;
				kv_push(aln_seed_t ,*p, *s_aln);
			}else { 
				push_chain(*av,aln_chain_t,*s_aln,aln_seed_t);
			}
		}
	}
}

/*
 *  Genrate chain :
 *        first step :  generate k-mer seed
 *        second step:  generate chains
 */

aln_chain_v *cm_mer2chain(const opt_t *opt , const seq_t *pseq)
{
	aln_chain_v  *av = malloc(sizeof(aln_chain_v));
	aln_seed_v   *kv_seed ;
	

	/* caculate  offset */
	int	offset =  cal_max_gap(opt,pseq->len) -  opt->l_seed ;
	
//      initialize struct aln_seed_v . Vector Size
	av->size = (int)((double)pseq->len/ opt->l_seed + 1.) ;
	round_size(av->size);
	aln_chain_v_init(*av);

	kv_seed = malloc(sizeof(aln_seed_v));
	kv_init(*kv_seed);
	av->offset = offset ;


	if(opt->verbose == 2){
		unit_sv_seed(kv_seed , opt);
		cm_mergechain(av,kv_seed);
		print_av_info(*av,opt);
		printf("************\n");
	}else  if(cm_meraln(kv_seed , pseq , opt)) {
		cm_mergechain(av,kv_seed);
		if( opt->verbose == 3){
			test_bns(*av ,pseq , opt);
		}
	}

	kv_destroy(*kv_seed);
	free(kv_seed);
	return  av;
}


/*
 * 
 * 	mark  flag which allow func cm_chain_extend to extend  the chain seed  
 *  	filter strategies is going to implement ..
 *
 */


void  mark_chain_se(aln_chain_v *av ,const opt_t *opt)
{
	int i , j ; 
	for( i = 0 ; i < av->n + 1 ; i++){
		aln_chain_t *at = av->a + i  ;
		//  calculate seedcov 
		at->seedcov = 0 ;
		for( j = 0  ;  j < at->n ; j++){
			aln_seed_t *p = at->a + j;
			at->seedcov += (p->query_e - p->query_b);
		}

		if(at->n >= 2 || at->seedcov >= 50 ) 
			at->flags =  ALLOW_EXTEND ;
		else	
			at->flags  = 0 ;
	}
}



#define  OVERLAP_REGION_SIZE	5
#define  SHORT_LENGTH		10
#define  MAX_LONG_LENGTH	50


/*
 *  @abstract:  if direct <= 0  , d is res's left side ... otherwise d is a's right side
 */
inline int  short_result( aln_seed_t *d ,  int direct , aln_res_t *res , const opt_t *opt )
{
	int ref_dist =  0 ;
	int read_dist = 0 ;
	int n_mis = 0 ;
	int n_ext = 0 ;

	d->flags &= DISALLOW_EXTEND ;

	if(direct){
		ref_dist =  d->ref_b  - res->ref_e ;
		read_dist =  d->query_b - res->query_e ;
	}else{
		ref_dist = res->ref_b -  d->ref_e ;
		read_dist = res->query_b - d->query_e ;
	}

	if(ref_dist == read_dist){
		n_mis = d->query_b - res->query_e ;
		res->n_mis += n_mis  ;
		res->app_score += (opt->a*d->len - opt->b*n_mis);
		//res->
	}else if ( abs(ref_dist) <  abs(read_dist)) {
		n_mis = ref_dist ? ref_dist : 0 ;
		res->n_mis += n_mis ;
		res->n_gap++ ;
		n_ext = (ref_dist < 0 ?  read_dist  -  ref_dist : read_dist) ;
		res->n_ext += n_ext ;
		res->app_score += (opt->a*d->len -  opt->b*n_mis - opt->o_del - opt->e_del*n_ext );
	}else {
		n_mis = read_dist ? read_dist : 0 ;
		res->n_mis +=  n_mis;
		res->n_gap++;
		n_ext = read_dist< 0 ?   ref_dist - read_dist: ref_dist;
		res->n_ext += n_ext ;
		res->app_score +=(opt->a*d->len - opt->b*n_mis - opt->o_ins - opt->e_ins*n_ext);
	}

	if(direct){
		res->ref_e = d->ref_e ;
		res->query_e =  d->query_e ;
	}else{
		res->ref_b  = d->ref_b ;
		res->query_b = d->query_b ; 
	}

	return  0 ;
}


/*
 * 	calculate linear penalty score , if penalty score > N  , this is a bad alignment.
 */

int  linear_aln( const opt_t *opt ,  const seq_t *seq , aln_seed_t *d , aln_res_t *res , int direct )
{
	int ref_dist =  0 ;
	int read_dist = 0 ;
	bwtint_t  ref_b  = 0 ; 
	int       query_b = 0 ;
	int  i  , k ;
	int64_t len ;

	unsigned char *bseq = NULL ;
	uint8_t *rseq = NULL ;
	kswr_t  r ;

	if(direct){
		ref_dist =  d->ref_b  - res->ref_e ;
		read_dist =  d->query_b - res->query_e ;
		ref_b =  res->ref_e ;
		query_b = res->query_e ;
	}else {
		ref_dist = res->ref_b -  d->ref_e ;
		read_dist = res->query_b - d->query_e ;
		ref_b = d->ref_e ;
		query_b = d->query_e ;
	}

	if(ref_dist == read_dist){

		char2nt4(&bseq,seq->len,seq->seq,0);
		seq_reverse(seq->len,bseq);
		rseq =  bns_get_seq(opt->fr->idx->bns->l_pac, opt->fr->idx->pac , ref_b , ref_b + ref_dist, &len );

		if(ref_dist != len )   printf("error\n");
		int  n_mis = 2 ;
		for( i = query_b + 1 ,  k = 1 ; i < query_b + read_dist - 1 && k < len -1   ; k++ , i++){
			if(bseq[i] !=  rseq[k])
				n_mis++;
		}
		if(opt->verbose == 5){

			for( i =  query_b   ;  i <  query_b + read_dist ; i++){
				printf("%c","ACGT"[bseq[i]]);
			}
			printf("\n");
			for( i = 0 ;  i  < len ; i++){
				printf("%c","ACGT"[rseq[i]]);

			}
			printf("\n");
			printf("%d\n",n_mis);


		}
		if((double)n_mis/read_dist > 0.5) {
			printf("exit\n");
			goto EXIT ;
		}else {
			res->n_mis += n_mis ;
			if(direct){
				res->query_e = d->query_e ;
				res->ref_e = d->ref_e ;

			}else{
				res->query_b = d->query_b ;
				res->ref_b =  d->ref_b ;
			}
			res->app_score += (opt->a*(d->len+read_dist) - opt->b*(n_mis));
			d->flags &= DISALLOW_EXTEND ;
		}
	}else {
		ref_b -= 25 ;
		query_b -= 25 ;
		if(ref_b < 0  || query_b < 0)  goto EXIT ;

		read_dist += 50;
		ref_dist += 50 ;
		if(ref_b + ref_dist >  (opt->fr->idx->bns->l_pac >> 1)  && query_b + read_dist > seq->len  )  goto EXIT;

		char2nt4(&bseq,seq->len,seq->seq,0);
		seq_reverse(seq->len,bseq);
		rseq =  bns_get_seq(opt->fr->idx->bns->l_pac, opt->fr->idx->pac , ref_b , ref_b + ref_dist, &len );
		if(opt->verbose == 6){

			for( i =  query_b   ;  i <  query_b + read_dist ; i++){
				printf("%c","ACGT"[bseq[i]]);
			}
			printf("\n");
			for( i = 0 ;  i  < len ; i++){
				printf("%c","ACGT"[rseq[i]]);

			}
			printf("\n");
		}
		r =  ksw_align2(read_dist , bseq + query_b , ref_dist , rseq , 5 ,opt->mat , opt->o_del ,opt->e_del ,opt->o_ins ,opt->e_ins , KSW_XSTART , 0);
		if(r.score >  35){
			d->flags &= DISALLOW_EXTEND ;
			if(direct){
				res->query_e = d->query_e ;
				res->ref_e = d->ref_e ;

			}else{
				res->query_b = d->query_b ;
				res->ref_b =  d->ref_b ;
			}
			res->app_score += (r.score - 50*opt->a) ; 
		}
	}
EXIT:
	free(bseq);
	free(rseq);
	return 0 ;

}
#define  MAX_BAND_TRY  2 
/*
 *  if  direct < 0 , left side  , else  right side. 
 */
inline int  extend_seed( const opt_t *opt , const  seq_t *seq ,  aln_res_t *res , int direct )
{
	uint64_t ref_b , ref_e ;
	uint8_t *rseq = NULL ;
	unsigned char *bseq = NULL ;
	int	query_b = 0  , query_e = 0 , i ;
	int64_t len ;
	int	max_off, aw = opt->w ;
	int	score =  -1  ;
	uint8_t  *qs , *rs ;
	if(direct){
		if(res->query_e == seq->len) 
			return  0 ;

		query_b = res->query_e;
		query_e = seq->len ;
		ref_b  =  res->ref_e ;
		ref_e  =  res->ref_e +  seq->len - res->query_e  >  opt->fr->idx->bns->l_pac << 1 ?  opt->fr->idx->bns->l_pac << 1 : res->ref_e + seq->len  - res->query_e ;
	}else {
	//  left side ..
		if(res->query_b == 0 )
			return 0 ;
		ref_b = res->ref_b - res->query_b  > 0  ? res->ref_b - res->query_b  :   0  ;
		ref_e =  res->ref_b ;
		query_b = 0 ;
		query_e = res->query_b ;
	}
	char2nt4(&bseq,seq->len,seq->seq,0);
	seq_reverse(seq->len,bseq);
	rseq =  bns_get_seq(opt->fr->idx->bns->l_pac, opt->fr->idx->pac , ref_b , ref_e , &len );
	if(direct){
		qs = bseq + query_b ;
		rs = rseq ;
	}else{
		qs = malloc((query_e-query_b)*sizeof(uint8_t));
		rs = malloc(len*sizeof(uint8_t));
		for( i = 0 ; i < query_e - query_b ; i++){
			qs[i] =  bseq[query_e - i -1];
		}
		for( i = 0 ;  i <  len ; i++){
			rs[i] = rseq[len-1-i]  ;
		}
	}
	if(opt->verbose == 7){
		if(direct) printf("right\n");
		else   printf("left\n");
		printf("%d , %d\n",res->query_b , res->query_e);

		for( i =  query_b   ;  i <  query_e  ; i++){
			printf("%c","ACGT"[bseq[i]]);
		}
		printf("\n");
		for( i = 0 ;  i  < len ; i++){
			printf("%c","ACGT"[rseq[i]]);

		}
		printf("\n");
		for( i = 0 ; i <  query_e - query_b ; i++){
			printf("%c","ACGT"[qs[i]]);
		}
		printf("\n");
		
		for( i = 0 ;  i  < len ; i++){
			printf("%c","ACGT"[rs[i]]);

		}
		printf("\n");
	}

	/*
	 *  maybe  direct ..
	 */

	int  qle, tle ,gtle ,gscore ;
	for( i = 0 ;  i  <  MAX_BAND_TRY  ; i++){
		int  pre =  score ;
		aw = aw << 1 ;
		score =  ksw_extend2(query_e-query_b , qs , len , rs , 5  , opt->mat , opt->o_del , opt->e_del ,opt->o_ins ,opt->e_ins , aw , opt->pen_clip5 , opt->zdrop ,res->app_score , &qle , &tle , &gtle , &gscore ,&max_off);
		if( score == pre  || max_off < (aw>>1)+(aw>>2)) break ;

	}
	printf("qle:%d\ttle:%d\tgtle:%d\tgscore:%dscore:%d\n",qle,tle,gtle,gscore,score);
	if(direct){
		if(gscore <= 0  ||  gscore <= score - opt->pen_clip5){
			res->query_e +=  qle ;
			res->ref_e +=  tle ;
		}else {
			res->query_e = seq->len;
			res->ref_e +=  gtle ;
		}

	}else {
		if(gscore <= 0  ||  gscore <= score - opt->pen_clip5){
			res->query_b -= qle ;
			res->ref_b -= tle ;
		}else {
			res->query_b = 0 ;
			res->ref_b -= gtle ;
		}
	}
	res->app_score = score ;

	if(opt->verbose == 7 ){
		print_res_info(res);
		printf("\n");
	}



	if(!direct){
		free(qs);
		free(rs);
	}

	free(rseq);
	free(bseq);
	return 0 ;
}

void cm_chain_extend(aln_chain_t *at, aln_res_v *rev ,const seq_t *seq , const opt_t *opt)
{
//	NOTE: p->ref_b   increasing ...	
	int i  , k  ;
// 
	for(i = 0 ; i < at->n ; i++){
		aln_seed_t *p = at->a + i ;
		p->flags  =  ALLOW_EXTEND ;
	}

	

	for( i = 0 ; i < at->n ; i++){
		aln_seed_t *p = at->a + i ;

		if(is_extend(*p)){
			//    push ... result ..
			aln_res_t   res ;
			//    reference  ID 
			res.ref_b  =  p->ref_b ;
			res.ref_e  =  p->ref_e ;
			res.query_b = p->query_b ;
			res.query_e = p->query_e ;
			res.app_score =  opt->a * ( p->query_e - p->query_b);
			res.n_mis = res.n_gap = res.n_ext =  res.sw_score = 0 ;
			//  left 
			for(  k = i-1   ;  k >= 0 ;  k--){
				aln_seed_t *res_left  = at->a + k ;
				if( res.ref_b + OVERLAP_REGION_SIZE <  res_left->ref_e  ||   res.query_b + OVERLAP_REGION_SIZE < res_left->query_b ){
				/* 
				 *  	When seeds have overlap region  , this seed can't be merge. 
				 *  	But we can't ignore these region  :  this reads maybe structure various
				 * */ 
					continue;
				}

				/*   when  long extend too long ,  these method's result  is not accurate. */

				if( res.ref_b - res_left->ref_e < OVERLAP_REGION_SIZE || res.query_b - res_left->query_e < OVERLAP_REGION_SIZE ){ 
					
					short_result(res_left,0,&res,opt);

				}else if ( res.ref_b - res_left->ref_e < MAX_LONG_LENGTH || res.query_b - res_left->query_e < MAX_LONG_LENGTH){
					
					linear_aln(opt , seq , res_left, &res , 0);
				}
			}
				// left side smith-waterman 
			extend_seed( opt , seq , &res , 0 );

			for(  k = i+1 ;  k < at->n ;  k++){

				aln_seed_t  *res_right = at->a + k ;
				if ( res_right->ref_b + OVERLAP_REGION_SIZE <  res.ref_e || res_right->query_b + OVERLAP_REGION_SIZE < res.query_e)  {

					continue;
				}
			

				if( res_right->ref_b - res.ref_e < OVERLAP_REGION_SIZE || res_right->query_b - res.query_e < OVERLAP_REGION_SIZE){
				/* 	
				 *  	short extend for reads left side. 
				 */
					short_result(res_right,1,&res,opt);
				} else if( res_right->ref_b - p->ref_e <  MAX_LONG_LENGTH|| res_right->query_b - p->query_e <  MAX_LONG_LENGTH){
				/* 	long  extend for reads left side. 
				 * 	too long  perform  smith-waterman ,
				 * 	perform  linear mapping  
				 */
					linear_aln(opt, seq , res_right , &res , 1);
					
				}

			}
			
			extend_seed( opt , seq , &res , 1 );
			kv_push( aln_res_t , *rev , res);
		}
	}
}

/*
 *        Those chain  is 'bad' , so  don't  need to extend . But we can't ignore every chain. 
 */

void  cm_disallow_extend(aln_chain_t *at , aln_res_v *rev , const opt_t *opt)
{
	int  i ;
	for( i = 0 ; i < at->n ; i++){
		aln_seed_t  *p = at->a + i ;
		
		aln_res_t  res  ;

		res.n_ext = 0 ;
		res.n_mis = 0 ;
		res.n_gap = 0 ;
		res.sw_score  =  0 ;
		res.rid =  0 ; 

		res.app_score = opt->a * (p->query_e - p->query_b );
		res.query_b =  p->query_b ;
		res.query_e =  p->query_e ;
		res.ref_b  = p->ref_b ;
		res.ref_e  = p->ref_e ;

		kv_push( aln_res_t , *rev , res);
	}
}

/*
 *	print sam stdout
 */

int  print_sam()
{
	return 0 ;
}

/*
 * 	chain ====> alignment result
 */

void  cm_chain2aln(aln_chain_v *av , aln_res_v *rev , const seq_t *seq ,const opt_t *opt ){
	int  i ;
	int  ret =  1 ; 
	//find_mismatch(seq->name,opt->l_seed,0);
	if(ret && opt->verbose == 4){
		printf("%s\n",seq->name);
	}
	for( i = 0 ; i < av->n + 1 ; i++){
		aln_chain_t *at = av->a + i  ;
		if(opt->verbose == 4 && ret){
			print_at_info(at);
		}
		if(is_extend(*at)) {
			cm_chain_extend(at,rev,seq,opt);
		}else	{
			cm_disallow_extend(at,rev,opt) ;
		}
	}
	if(opt->verbose == 4 && ret){
		print_resv_info(rev);
	}

}

int  cm_chain_core(const opt_t *opt , const mseq_t *mseq , aln_res_v  *rev)
{
	int  i ; 

	for( i = 0 ; i  <  mseq->n_seq ; i++){
		
		seq_t *p =  mseq->seq + i ;
		aln_chain_v *av = cm_mer2chain(opt,p);

		if(opt->verbose == 1 )
			test_pos(p->name , *av , 0 , opt->l_seed , opt);

/*
 *		mark chains  which not need verification.
 *
 */
		mark_chain_se(av,opt);

/*
 *		chain extend by linaer alignment and  smith-waterman extension.
 */
		cm_chain2aln(av,&rev[i],p,opt);


		if(opt->verbose == 5)  //  extend test 
			;  

		
/*
 * 	step4: Fmeas filter ...
 */

		aln_chain_v_free(*av);
		free(av);
	}

	return 0 ;
}



int   multi_seq_get(seqdb_t db, mseq_t *s){

	int i = 0 ;
	s->n_seq = 0; 
	s->seq = malloc(s->n_need*sizeof(seq_t));
	for( i = 0 ; i < s->n_need ; i++){
		seq_t *p;
		seqdb_get(db,&p);
		if(!p) break;
		s->seq[i] = *p ;
		free(p);
	}
	s->n_seq = i ;
	s->total_num += i ;
	return 0;
}

void  multi_seq_release(seq_t *seq){

	if(seq){
		free(seq->seq);
		free(seq->name);
		free(seq->add);
		free(seq->qual);
	}
}

void  multi_seq_free(mseq_t *s, int sel)
{
	int  i , j ;
	for( i  = 0 ;  i  < sel ; i++ ){
		mseq_t *p =  s + i ;
		for( j = 0 ; j  < p->n_seq ; j++){
			seq_t  *seq =  p->seq + j ;
			multi_seq_release(seq) ;

		}
		free(p->seq);
	}
}

/*
 * 	this section select good alignment , then  output.
 */

int	cm_se( const opt_t *opt , aln_res_v *rev)
{
	print_sam();

	return 0 ;
}


int  aln_core(const opt_t *opt)
{
	mseq_t  *mseq;
	int  i ;
//     store alignment result 	
	aln_res_v  *rev[2] ;
	
	if(opt->flag & F_PE){
		mseq = calloc(2,sizeof(mseq_t));
		mseq[0].n_need = opt->n_need;
		mseq[1].n_need = opt->n_need;
		rev[0]  = malloc(opt->n_need*sizeof(aln_res_v));
		rev[1]  = malloc(opt->n_need*sizeof(aln_res_v));
		for( i = 0 ; i < opt->n_need ; i++){
			kv_init(rev[0][i]);
			kv_init(rev[1][i]);
		}

	}else  {
		mseq =  calloc(1,sizeof(mseq_t));
		mseq[0].n_need = opt->n_need ;
		rev[0]  = malloc(opt->n_need*sizeof(aln_res_v));
		for( i = 0 ; i < opt->n_need ; i++)
			kv_init(rev[0][i]);
		
	}

	do{

		if(opt->flag & F_PE){
			multi_seq_get(opt->fs[0].seqdb,&mseq[0]);
			multi_seq_get(opt->fs[1].seqdb,&mseq[1]);
			if(mseq[0].n_seq != mseq[1].n_seq){
				fprintf(stderr,"two sequece file is not pair file\n");
				return EXIT_FAILURE ;
			}
		/*
		 *  check reads name..  but DNAA-cigar is wrong
		 */
		}else  multi_seq_get(opt->fs->seqdb,mseq);
		
		
		if( mseq[0].n_seq == 0 ) {
			multi_seq_free(mseq,1 + (opt->flag&F_PE));
			break;
		}

		if(opt->flag & F_PE){
			cm_chain_core(opt,&mseq[0],rev[0]);
			cm_chain_core(opt,&mseq[1],rev[1]);
		}else   cm_chain_core(opt,&mseq[0],rev[0]);

		if(opt->flag & F_PE)  
			cm_pe(opt,rev);
		else	
			cm_se(opt,rev[0]);
		if(opt->flag & F_PE){
			for( i = 0 ; i < opt->n_need ;i++){
				rev[0][i].n = rev[1][i].n = 0 ;
			}
		}else{
			for( i = 0 ; i < opt->n_need ;i++){
				rev[0][i].n = 0 ;
			}

		}
	
		
		multi_seq_free(mseq, 1 + (opt->flag&F_PE));
		fprintf(stderr,"%ld reads have been processed\n",mseq[0].total_num);

	}while(1);

	for( i = 0 ; i < opt->n_need; i++){
		kv_destroy(rev[0][i]);
	}

	if((opt->flag&F_PE)){
		for( i = 0 ; i < opt->n_need ; i++)
			kv_destroy(rev[1][i]);
	}

	free(rev[0]);
	free(rev[1]);
	free(mseq);
	return EXIT_SUCCESS;
}
