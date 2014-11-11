#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "aln.h"
#include "kvec.h"



void  print_av_info(aln_chain_v av , const opt_t *opt){
	int j , k ;
	for( j = 0 ;  j  <  av.n + 1 ; j++){
		aln_chain_t *p  = av.a + j ;
		for( k = 0 ; k < p->n ; k++ ){
			aln_seed_t *sp =   p->a + k ;
			int is_rev;
			bwtint_t pos = bns_depos(opt->fr->idx->bns, sp->ref_b > opt->fr->idx->bns->l_pac ? sp->ref_e - 1 : sp->ref_b, &is_rev);  
			printf("%ld\t%ld\t",pos,pos+sp->ref_e-sp->ref_b);
		}
	}
	printf("\n");
}

void  print_at_info(aln_chain_t *at ){

	int  i ;
	for(  i  = 0 ; i  <  at-> n ; i++){
		aln_seed_t  *p  =  at->a +  i ;
		printf("# qb:%d qe:%d tb:%ld te:%ld\n",p->query_b , p->query_e , p->ref_b , p->ref_e);

	}
	printf("#\n");

}

int  test_bns(const aln_chain_v av , const seq_t *seq ,const opt_t *opt )
{
	int i ,j , k , l;
	for( i = 0 ; i < av.n + 1 ; i++){
		aln_chain_t *p =  av.a + i ;
		for( j = 0 ; j < p->n ; j++){
			aln_seed_t *q = p->a + j ;
			unsigned char *bseq ;
			int64_t  rlen ;
			uint8_t	 *rseq = NULL ;  //  reference sequence ..
			char2nt4(&bseq,seq->len,seq->seq,0);
			seq_reverse(seq->len,bseq);
			rseq = bns_get_seq(opt->fr->idx->bns->l_pac , opt->fr->idx->pac , q->ref_b , q->ref_e , &rlen);
			for( k = seq->len - q->query_e  ; k < seq->len - q->query_b ; k++){
				printf("%c","ACGT"[bseq[k]]);
			}
			printf("\n");
			for ( l = 0 ; l < rlen ; l++)
				printf("%c","ACGT"[rseq[l]]);
			printf("\n");
			free(bseq);
			free(rseq);
		}
	}
	return 0 ;
}

int  key_cmp(const void *p1  , const void *p2)
{
	const  key_pos_t  *ss1 = p1 ;
	const  key_pos_t  *ss2 = p2 ;
	return  (ss2->beg > ss1->beg) - (ss1->beg > ss2->beg);
}

typedef struct {
	int  pos[2];
	int  strand[2];
	key_pos_v  *kpv;
} DNAA_info ;

DNAA_info *gain_test_name_info(char *name , int l_seed , int sel)
{
	char  *tmp =  strstr(name,"rand");
	DNAA_info   *info  = NULL  ; 
	int   i , n_cigar ;

	if(tmp)  return info ;

	info = malloc(sizeof(DNAA_info));

	tmp = name ;

	for( i = 0 ; i < 7 ; i++){
		tmp =  strstr(tmp,"_");
		tmp++;
		switch(i){
			case 0:  info->pos[0] = atoi(tmp); break;
			case 1:  info->pos[1] = atoi(tmp); break;
			case 2:  info->strand[0] = atoi(tmp); break;
			case 3:  info->strand[1] = atoi(tmp); break;

		}
		
	}
//	printf("1err:%d\t",atoi(tmp)); 
	tmp =  strstr(tmp,":");
	tmp++;
//	printf("1sub:%d\t",atoi(tmp));
	tmp =  strstr(tmp,":");
	tmp++;
//	printf("indel:%d\t",atoi(tmp));
	
	tmp =  strstr(tmp,"_");
	tmp++;
	
//	printf("2err:%d\t",atoi(tmp)); 
	tmp =  strstr(tmp,":");
	tmp++;
//	printf("2sub:%d\t",atoi(tmp));
	tmp =  strstr(tmp,":");
	tmp++;
//	printf("indel:%d\t",atoi(tmp));
	
	tmp =  strstr(tmp,"_");
	tmp++;
	
	tmp =  strstr(tmp,"_");
	tmp++;

	char *p ;
	n_cigar = 0 ;

	info->pos[0]--;
	info->pos[1]--;

	for( p = tmp ; *p != '/' ; p++)
		if(!isdigit(*p)) n_cigar++;

	cigar_t *cigar =  malloc(sizeof(cigar_t)*n_cigar);
	for( i = 0 ;  i < n_cigar ; i++,tmp++){
		cigar[i].len = strtol(tmp,&tmp,10);
		cigar[i].op =  *tmp;
	//	printf("%d%c",cigar[i].len,*tmp);
	}
	info->kpv = malloc(sizeof(key_pos_v));
	info->kpv->val = calloc(n_cigar,sizeof(key_pos_t));
	info->kpv->n  = 0 ;

	int	pos =  info->pos[sel];
	
	for( i = 0 ; i < n_cigar ; i++){
		cigar_t *c = cigar + i ;
		switch(c->op){
			case 'M': 
				if(c->len >= l_seed) {
				//	printf("%d\t%d\t",pos[sel],pos[sel]+c->len); 
					key_pos_t  *val =  info->kpv->val + info->kpv->n ;
					val->beg =  pos ;
					val->end =  pos + c->len ;
					info->kpv->n++ ;
				}
				pos += c->len ;
				break;
			case 'U':  pos+= c->len ; break ;
			case 'I':  break ; 
			case 'D':  pos+= c->len ; break ;
		}

	}

	free(cigar);

	return info ;
}

/*
 * 	unit for chain's pos , overlap exmaple 21M1I70M   return  22M  70M .
 */

int	test_pos(char *name ,  const aln_chain_v  av , int sel , int l_seed ,const opt_t *opt)
{
	DNAA_info *info = gain_test_name_info(name,l_seed,sel);
	if(!info) return 0 ;
	if(!info->strand[sel]) qsort(info->kpv->val , info->kpv->n , sizeof(key_pos_t), key_cmp );
	int i , j , k ;
	for( j = 0 ;  j  <  av.n + 1 ; j++){
		aln_chain_t *p  = av.a + j ;
		for( k = 0 ; k < p->n ; k++ ){
			aln_seed_t *sp =   p->a + k ;
			int is_rev;
			bwtint_t pos = bns_depos(opt->fr->idx->bns, sp->ref_b > opt->fr->idx->bns->l_pac ? sp->ref_e - 1 : sp->ref_b, &is_rev);  
		//	printf("%ld\t%ld\t",pos,pos+sp->ref_e-sp->ref_b);
			for( i = 0 ; i < info->kpv->n  ; i++){
				key_pos_t  *val =  info->kpv->val + i ;
				if(val->beg == pos &&  pos+sp->ref_e-sp->ref_b == val->end){
					val->right = 1 ;
				}
			}
		}
	}
	int  err = 0 ;
	for( i = 0 ; i < info->kpv->n  ; i++){
		key_pos_t  *val =  info->kpv->val + i ;
		if(val->right == 0)  err = 1 ;
	}
	if(err){
		for( i = 0 ; i < info->kpv->n  ; i++){
			key_pos_t  *val =  info->kpv->val + i ;
			printf("%d\t%d\t",val->beg,val->end);
		}
		for( j = 0 ;  j  <  av.n + 1 ; j++){
			aln_chain_t *p  = av.a + j ;
			for( k = 0 ; k < p->n ; k++ ){
				aln_seed_t *sp =   p->a + k ;
				int is_rev;
				bwtint_t pos = bns_depos(opt->fr->idx->bns, sp->ref_b > opt->fr->idx->bns->l_pac ? sp->ref_e - 1 : sp->ref_b, &is_rev);  
				printf("%ld\t%ld\t",pos,pos+sp->ref_e-sp->ref_b);
			}
		}
		printf("%s",name);
	}
	free(info->kpv->val);
	free(info->kpv);
	free(info);
	return 0 ;
}


/*
 * 	unit for the repeat same coordinate on difference seeds
 * 	generate  this seed.  verbose == 2 
 */

void  unit_sv_seed(aln_seed_v *kv_seed , const opt_t *opt)
{
	int  i ;
	aln_seed_t  s_aln ;
	for( i = 0 ; i < 10 ; i++){
		s_aln.ref_b = 100 - i ;
		s_aln.ref_e = 100 + opt->l_seed - i ;
		s_aln.query_b = i ;
		s_aln.query_e = i+opt->l_seed ;
		kv_push(aln_seed_t , *kv_seed ,s_aln);
	#if 1
		s_aln.ref_b = 100 - i ;
		s_aln.ref_e = 100 + opt->l_seed - i ;
		s_aln.query_b = i + 50 ;
		s_aln.query_e = i +opt->l_seed + 50;
		kv_push(aln_seed_t , *kv_seed ,s_aln);
	#endif 
	}
	#if 1
	i += opt->l_seed;
	for( ; i < 50 ; i++){
		s_aln.ref_b = 100 - i ;
		s_aln.ref_e = 100 + opt->l_seed - i ;
		s_aln.query_b = i ;
		s_aln.query_e = i+opt->l_seed ;
		kv_push(aln_seed_t , *kv_seed ,s_aln);
	}
	#endif
	for( i = 0 ; i < kv_seed->n ; i++){
		aln_seed_t  *p = kv_seed->a + i ;
		printf("%ld\t%ld\t%d\t%d\n",p->ref_b , p->ref_e , p->query_b , p->query_e);

	}
	printf("****************************************\n");


}
/*
 * 	unit for two seed short merge ..
 */
int	unit_extend_1( char *name , aln_res_v  *rev , int sel , int l_seed , const opt_t *opt)
{
	DNAA_info *info = gain_test_name_info(name,l_seed,sel);
	if(!info) return 0 ;

	int  i ;
	for( i = 0 ; i < rev->n ; i++){
		aln_res_t  *p =  rev->a + i ; 
		printf("%d\n",p->app_score);
	}

	free(info->kpv->val);
	free(info->kpv);
	free(info);
	return 0 ;
}
