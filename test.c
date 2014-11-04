#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "aln.h"



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

void  print_at_info(aln_chain_t *at , const opt_t *opt){

	int  i ;
	for(  i  = 0 ; i  <  at-> n ; i++){
		aln_seed_t  *p  =  at->a +  i ;
		printf("# qb:%d qe:%d tb:%ld te:%ld",p->query_b , p->query_e , p->ref_b , p->ref_e);

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


int	test_pos(char *name ,  const aln_chain_v  av , int sel , int l_seed ,const opt_t *opt)
{
	char  *tmp =  strstr(name,"rand");
	int   i , n_cigar ;
	int   pos[2] , strand[2] ;

	if(tmp){
		return 0 ;
	}
	tmp = name ;

	for( i = 0 ; i < 7 ; i++){
		tmp =  strstr(tmp,"_");
		tmp++;
		switch(i){
			case 0:  pos[0] = atoi(tmp); break;
			case 1:  pos[1] = atoi(tmp); break;
			case 2:  strand[0] = atoi(tmp); break;
			case 3:  strand[1] = atoi(tmp); break;

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

	pos[0]--;
	pos[1]--;

	for( p = tmp ; *p != '/' ; p++)
		if(!isdigit(*p)) n_cigar++;

	cigar_t *cigar =  malloc(sizeof(cigar_t)*n_cigar);
	for( i = 0 ;  i < n_cigar ; i++,tmp++){
		cigar[i].len = strtol(tmp,&tmp,10);
		cigar[i].op =  *tmp;
	//	printf("%d%c",cigar[i].len,*tmp);
	}
//	printf("\t");
	key_pos_v  *kpv = malloc(sizeof(key_pos_v));
	kpv->val = calloc(n_cigar,sizeof(key_pos_t));
	kpv->n  = 0 ;

	
	for( i = 0 ; i < n_cigar ; i++){
		cigar_t *c = cigar + i ;
		switch(c->op){
			case 'M': 
				if(c->len >= l_seed) {
				//	printf("%d\t%d\t",pos[sel],pos[sel]+c->len); 
					key_pos_t  *val =  kpv->val + kpv->n ;
					val->beg =  pos[sel] ;
					val->end =  pos[sel] + c->len ;
					kpv->n++ ;
				}
				pos[sel]+= c->len ;
				break;
			case 'U':  pos[sel]+= c->len ; break ;
			case 'I':  break ; 
			case 'D':  pos[sel]+= c->len ; break ;
		}

	}
	if(!strand[sel]) qsort(kpv->val , kpv->n , sizeof(key_pos_t), key_cmp );
	/* 		
	}
	*/ 
	int j , k ;
	for( j = 0 ;  j  <  av.n + 1 ; j++){
		aln_chain_t *p  = av.a + j ;
		for( k = 0 ; k < p->n ; k++ ){
			aln_seed_t *sp =   p->a + k ;
			int is_rev;
			bwtint_t pos = bns_depos(opt->fr->idx->bns, sp->ref_b > opt->fr->idx->bns->l_pac ? sp->ref_e - 1 : sp->ref_b, &is_rev);  
		//	printf("%ld\t%ld\t",pos,pos+sp->ref_e-sp->ref_b);
			for( i = 0 ; i < kpv->n  ; i++){
				key_pos_t  *val =  kpv->val + i ;
				if(val->beg == pos &&  pos+sp->ref_e-sp->ref_b == val->end){
					val->right = 1 ;
				}
			}
		}
	}

	int  err = 0 ;
	for( i = 0 ; i < kpv->n  ; i++){
		key_pos_t  *val =  kpv->val + i ;
		if(val->right == 0)  err = 1 ;
	}
	if(err){
		for( i = 0 ; i < kpv->n  ; i++){
			key_pos_t  *val =  kpv->val + i ;
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


	free(cigar);
	free(kpv->val);
	free(kpv);


	return 0 ;
}
