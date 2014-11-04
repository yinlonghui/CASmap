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


#define flt_fuc(a,b) ((a).ref_b < (b).ref_b)
KSORT_INIT(cm_seed_flt,aln_seed_t,flt_fuc)



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
		if(bwt_match_exact_cm( opt->fr->idx->bwt , i , opt->l_seed ,bseq  , &k , &l)){ // set  MAX_INTERVAL  , to filter redundancy k-mer hit..
			bwtint_t  j ,pos;
			for( j = k ;  j < l + 1 ; j++){
				pos = bwt_sa(opt->fr->idx->bwt,j);
				s_aln.ref_b =  pos , s_aln.ref_e = pos + opt->l_seed ;
				s_aln.query_b = i  , s_aln.query_e =  i + opt->l_seed ;
	//			printf("query_b %d query_e %d ref_b %ld ref_e %ld \n",s_aln.query_b,s_aln.query_e,s_aln.ref_b,s_aln.ref_e);
				kv_push(aln_seed_t,*kv_seed,s_aln);
			}
		}
	}
	free(bseq);
	return  kv_size(*kv_seed);
}

void   cm_mergechain(aln_chain_v *av ,  aln_seed_v  *kv_seed)
{
	int  j ; 
	
	aln_seed_t  s_aln;
//       first , regular merge seed .
	ks_introsort(cm_seed_flt,kv_size(*kv_seed),kv_seed->a);
	s_aln = kv_A(*kv_seed,0)	;
	chain_at_add(av->a[av->n],s_aln);
	for ( j = 1 ;  j < kv_size(*kv_seed) ; j++){
		aln_seed_t p_aln = kv_A(*kv_seed,j);
		s_aln =  chain_av_last(*av);
		int l_ref = s_aln.ref_e -  p_aln.ref_e ;
		int l_read = p_aln.query_b  - s_aln.query_b ;
		if( l_read == l_ref && s_aln.ref_e >= p_aln.ref_b ){
			s_aln.ref_e  =  p_aln.ref_e ;
			s_aln.query_b = p_aln.query_b  ;
			chain_av_last_c(*av,s_aln);
		}else if( p_aln.ref_e +  av->offset >=  s_aln.ref_b){
			chain_at_add(av->a[av->n],p_aln);
		}else {
			chain_av_add(*av,p_aln);
		}

	}
/*
*	  check whether chain have unmerge seed.
*	  eg:  ref: ------ACCCTGATCCCTGA------
*	  	reads     ACCCTGATCCCTGACCCTGATCCC
*	 example:   tandemCNV 
*	 test unite ...
*/ 

//     second , sv seed merge 
}


aln_chain_v *cm_mer2chain(const opt_t *opt , const seq_t *pseq)
{
	aln_chain_v  *av = malloc(sizeof(aln_chain_v));
	aln_seed_v   *kv_seed ;
	

	int	offset = pseq->len + cal_max_gap(opt,pseq->len) -  opt->l_seed ;
	
//      initialize struct aln_seed_v 	
	av->size = (int)((double)pseq->len/ opt->l_seed + 1.) ;
	round_size(av->size);
	aln_chain_v_init(*av);

	kv_seed = malloc(sizeof(aln_seed_v));
	kv_init(*kv_seed);
	av->offset = offset ;
	if(cm_meraln(kv_seed , pseq , opt))   cm_mergechain(av,kv_seed);
	kv_destroy(*kv_seed);
	free(kv_seed);
	return  av;
}




void  mark_chain_se(aln_chain_v *av)
{
	int i  ; 
	for( i = 0 ; i < av->n ; i++){
		aln_chain_t *at = av->a + i  ;
		at->flags =  ALLOW_EXTEND ;
		//  before extend , calculate score ..  
		//  filter strategies is going to implement ..
	}
}

void cm_chain_extend(aln_chain_t *at  ,const opt_t *opt)
{
	int  i = 0 ;
	for( i = 0 ; i < at->n ; i++){
		// short ..
		
		// long
		
		// smith-waterman 


	}


}

void  cm_disallow_extend(aln_chain_t *at , const opt_t *opt)
{
	int  i ;
	for( i = 0 ; i < at->n ; i++){
	//	res.sw_score  = 0 ;
		;
	}
}

int  print_sam()
{
	return 0 ;
}


void  cm_chain2aln(aln_chain_v *av , const opt_t *opt ){
	int  i ;
	for( i = 0 ; i < av->n ; i++){
		aln_chain_t *at = av->a + i  ;
		if(is_extend(at)) cm_chain_extend(at,opt);
		else	cm_disallow_extend(at,opt) ;
	}
}

int  cm_chain_core(const opt_t *opt , const mseq_t *mseq , aln_res_v  *rev)
{
	int  i ; 

	for( i = 0 ; i  <  mseq->n_seq ; i++){
		
		seq_t *p =  mseq->seq + i ;
		aln_chain_v *av = cm_mer2chain(opt,p);
/*
 *		mark chains  which not need verification.
 *
 */
//		mark_chain_se(&av);
		if(opt->verbose == 4 )test_pos(p->name , *av , 0 , opt->l_seed , opt);

/*
 *		chain extend by linaer alignment and  smith-waterman extension.
 */
//		cm_chain2aln(&av,opt);
		
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
	aln_res_v  *rev ;
	
	if(opt->flag & F_PE){
		mseq = calloc(2,sizeof(mseq_t));
		mseq[0].n_need = opt->n_need;
		mseq[1].n_need = opt->n_need;
		rev  = malloc(2*sizeof(aln_res_v));
		kv_init(rev[0]);
		kv_init(rev[1]);

	}else  {
		mseq =  calloc(1,sizeof(mseq_t));
		mseq[0].n_need = opt->n_need ;
		rev  = malloc(sizeof(aln_res_v));
		kv_init(*rev);
	}

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
			cm_chain_core(opt,&mseq[0],&rev[0]);
			cm_chain_core(opt,&mseq[1],&rev[1]);
		}else   cm_chain_core(opt,&mseq[0],rev);

		if(opt->flag & F_PE)  cm_pe(opt,rev);
		else	cm_se(opt,rev);

		multi_seq_free(mseq, 1 + (opt->flag&F_PE));
		fprintf(stderr,"%ld reads have been processed\n",mseq[0].total_num);

	}while(1);
	for( i = 0 ; i < 1 + (opt->flag&F_PE) ; i++){
		kv_destroy(rev[i]);
	}
	free(rev);
	free(mseq);

	return EXIT_SUCCESS;
}
