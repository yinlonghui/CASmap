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



static inline int cal_max_gap (const opt_t *opt , int len)
{	
	int  l_del = (int)((double)(len*opt->a - opt->o_del)/opt->e_del + 1.);
	int  l_ins = (int)((double)(len*opt->a - opt->o_ins)/opt->e_ins + 1.);
	int  l =  l_del  >  l_ins  ?  l_del : l_ins ;
	l = l > 1 ?  l : 1 ;
	return  l < opt->w << 1 ? l : opt->w << 1;
}

/*
 * 	This alignment mode allow mismatch so that we must using the stack.
 */

#define flt_fuc(a,b) ((a).ref_b < (b).ref_b)
KSORT_INIT(cm_seed_flt,aln_seed_t,flt_fuc)
aln_chain_v cm_inser_seed(const opt_t *opt , const seq_t *pseq)
{

//	int len = 0 ,beg = 0 ;
	
	unsigned char *bseq ;

	int	i  =  0  ;


	kvec_t(aln_seed_t) kv_seed;
	aln_seed_t  s_aln;

	aln_chain_v  av;
	
	kv_init(kv_seed);

	av.size = (int)((double)pseq->len/ opt->l_seed + 1.) ;
	round_size(av.size);
	av.offset =  pseq->len + cal_max_gap(opt,pseq->len) -  opt->l_seed ;
	aln_chain_v_init(av);

/*
 *	step1:	Generate seed  insert to chain 
 */

//	len =  opt->l_seed ;
	char2nt4(&bseq,pseq->len,pseq->seq,0); seq_complement(pseq->len ,bseq);
//	exact match  int bwt_match_exact_cm(const bwt_t *bwt,int begin , int len , const unsigned char *str, bwtint_t *pk , bwtint_t *pl)
	for( i  = 0 ;  i  <  pseq->len - opt->l_seed + 1 ; i++){
		bwtint_t  k , l;
		k = 0 , l = opt->fr->idx->bwt->seq_len ;
		if(bwt_match_exact_cm( opt->fr->idx->bwt , i , opt->l_seed ,bseq  , &k , &l)){ // set  MAX_INTERVAL  , to filter redundancy k-mer hit..
			bwtint_t  j ,pos;
			for( j = k ;  j < l + 1 ; j++){
				pos = bwt_sa(opt->fr->idx->bwt,j);
				s_aln.ref_b =  pos , s_aln.ref_e = pos + opt->l_seed ;
				s_aln.query_b = i  , s_aln.query_e =  i + opt->l_seed ;
	//			printf("query_b %d query_e %d ref_b %ld ref_e %ld \n",s_aln.query_b,s_aln.query_e,s_aln.ref_b,s_aln.ref_e);
				kv_push(aln_seed_t,kv_seed,s_aln);
			}
		}
	}

	int j = 0  ;
/*
 *      step2:  sorting  the seed ..
 *      There are many algorithm to choose..
 */
	ks_introsort(cm_seed_flt,kv_size(kv_seed),kv_seed.a);
	/*
	 *  merge to one seed
	 */ 
	if(!kv_size(kv_seed)){
		kv_destroy(kv_seed);
		free(bseq);
		return av;
	}
	
	s_aln = kv_A(kv_seed,0)	;
	chain_at_add(av.a[av.n],s_aln);
/*
 * 	merge  k-mer seed to  longer seed 
 */
	
	for ( j = 1 ;  j < kv_size(kv_seed) ; j++){
		aln_seed_t p_aln = kv_A(kv_seed,j);
		s_aln =  chain_av_last(av);
		int l_ref = s_aln.ref_e -  p_aln.ref_e ;
		int l_read = p_aln.query_b  - s_aln.query_b ;
		if( l_read == l_ref && s_aln.ref_e >= p_aln.ref_b ){
			s_aln.ref_e  =  p_aln.ref_e ;
			s_aln.query_b = p_aln.query_b  ;
			chain_av_last_c(av,s_aln);
		}else if( p_aln.ref_e +  av.offset >=  s_aln.ref_b){
			chain_at_add(av.a[av.n],p_aln);
		}else {
			chain_av_add(av,p_aln);
		}

	}
/*
*	  check whether chain have unmerge seed.
*	  eg:  ref: ------ACCCTGATCCCTGA------
*	  	reads     ACCCTGATCCCTGACCCTGATCCC
*	 example:   tandemCNV 
*	 test unite ...
*/ 

	kv_destroy(kv_seed);
	free(bseq);
	return  av;
}


/*
 * 	smith -waterman  alignment struct ...
 */
typedef struct {
	int query_b , query_e ;
	bwtint_t ref_b ,  ref_e ;

} aln_sw_t ;

typedef  struct {
	int m , n ;
	aln_sw_t   *a ;
} aln_sw_v ;

void print_swv(aln_sw_v  asw)
{
	int i ;
	for( i = 0 ; i < asw.n ; i++){
		aln_sw_t  *p  =  asw.a + i ;
		printf("rb = %ld  re = %ld  qb = %d  qe = %d \n",p->ref_b , p-> ref_e , p->query_b ,p->query_e);
	}

}
/*
 *	step3: Chains can extend by Smith-waterman , The smith-waterman function in ksw.c can be reused.
 */

#define MAX_BAND_TRY 2 
#define CM_SHORT_LEN 200
#define CM_SHORT_EX  50

/*
 *   The short alignment based on sse2  don't  increase the accuracy. 
 */

void  mark_chain_se(aln_chain_v *av)
{
	int i  ; 
	for( i = 0 ; i < av->n ; i++){
		aln_chain_t *at = av->a + i  ;
		at->extend =  1 ;
		//  before extend , calculate score ..  
		//  filter strategies is going to implement ..
	}
}


void  approximate_align(aln_chain_t *at)
{
	int i;
	for( i = 0 ; i < at->n ; i++){



	}
}


void  extend_align(aln_chain_t *at)
{
	;
}
void  chain2result(aln_chain_t *at)
{
	;
}

int  print_sam()
{
	return 0 ;
}

#define  cal_score(at)  (at->score =  0) 

void  chain_extend(aln_chain_v *av){
	int  i ;
	for( i = 0 ; i < av->n ; i++){
		aln_chain_t *at = av->a + i  ;
		if(at->extend){
			//  approximately linear alignment 
			approximate_align(at);	
			//  extend by using smith-waterman 
			extend_align(at);
			//  after extend , calculate score ..
			cal_score(at);
		}else   chain2result(at);

	}

}

int  cm_chain_core(const opt_t *opt , const mseq_t *mseq)
{
	int  i ; 
	for( i = 0 ; i  <  mseq->n_seq ; i++){
		
		seq_t *p =  mseq->seq + i ;
		aln_chain_v av = cm_inser_seed(opt,p);
		aln_sw_v    asw;
		kv_init(asw);

		if(opt->verbose == 2 ){
			printf("%s\t",p->name);
			print_av_info(av,opt);
		}
		if(opt->verbose == 3){
			printf("%s\n",p->name);
			test_bns(av,p,opt);
		}
		if(opt->verbose == 4){
			test_pos(p->name ,av ,0 ,opt->l_seed ,opt);
		}

/*
 *		mark chains  which not need verification.
 *
 */
		mark_chain_se(&av);

/*
 *		chain extend by linaer alignment and  smith-waterman extension.
 */
		chain_extend(&av);


		
		if(opt->verbose == 6){
			if(asw.n)  printf("%s",p->name);
			print_swv(asw);
		}
		
/*
 * 	step4: Fmeas filter ...
 */
		kv_destroy(asw);
		aln_chain_v_free(av);
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

int	cm_se()
{
	print_sam();

	return 0 ;
}


int  aln_core(opt_t *opt)
{
	mseq_t  *mseq;
	
	if(opt->flag & F_PE){
		mseq = calloc(2,sizeof(mseq_t));
		mseq[0].n_need = opt->n_need;
		mseq[1].n_need = opt->n_need;
	}else  {
		mseq =  calloc(1,sizeof(mseq_t));
		mseq[0].n_need = opt->n_need ;
	}

/*
 *	The optimizing strategy is establishing the cache to access 
 */

	do{

		if(opt->flag & F_PE){
			multi_seq_get(opt->fs[0].seqdb,&mseq[0]);
			multi_seq_get(opt->fs[1].seqdb,&mseq[1]);
			if(mseq[0].n_seq != mseq[1].n_seq){
				fprintf(stderr,"two sequece file is not pair file\n");
				return EXIT_FAILURE ;
			}
		}else  multi_seq_get(opt->fs->seqdb,mseq);
		
		
		if( mseq[0].n_seq == 0 ) {
			multi_seq_free(mseq,1 + (opt->flag&F_PE));
			break;
		}

		if(opt->flag & F_PE){
			cm_chain_core(opt,&mseq[0]);
			cm_chain_core(opt,&mseq[1]);
		}else   cm_chain_core(opt,&mseq[0]);

		if(opt->flag & F_PE)  cm_pe();
		else	cm_se();

		multi_seq_free(mseq, 1 + (opt->flag&F_PE));
		fprintf(stderr,"%ld reads have been processed\n",mseq[0].total_num);

	}while(1);


	free(mseq);
	return EXIT_SUCCESS;
}
