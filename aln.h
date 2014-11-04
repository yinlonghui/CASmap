#ifndef __ALN_H
#define __ALN_H
#include "seq.h"
#include "bwa.h"
#define  F_PE  1

typedef struct {
	seqdb_t	seqdb;
	char	*fn;
}  file_seq_t;

typedef struct{
	char	*fn;
	bwaidx_t *idx;
}  file_ref_t;

typedef struct {
	file_seq_t *fs;
	file_ref_t *fr;
	int8_t mat[25];
	int  a , b;
	int  o_del , o_ins ;
	int  e_ins , e_del ;
	int  pen_clip5 , pen_clip3 ;
	int  zdrop ;
	int  w;
	int  n_need ;
	int  l_seed ;
	int  CC_stack;
	int  verbose;
	int  flag;
} opt_t ;

typedef struct {
	seq_t	*seq;
	uint64_t  n_need;
	uint64_t  total_num ;
	uint64_t  n_seq;
} mseq_t ;

typedef struct {
	bwtint_t ref_b,   ref_e;
	int	 query_b, query_e;
	int	 flags ;
} aln_seed_t;

#define  ALLOW_EXTEND  1 
#define  HAVE_OVERLAP  2

#define  is_extend(at) (at->flags & ALLOW_EXTEND)
#define  is_overlap(at) (at->flags & HAVE_OVERLAP)

typedef struct {
	aln_seed_t  *a ;
	int	flags;
	int	m,n;
} aln_chain_t ;

typedef struct {
	aln_chain_t *a;
	int     m,n,size,offset;
} aln_chain_v;

typedef struct {
	aln_seed_t *a  ;
	int n , m ;
} aln_seed_v ;

typedef struct{
	bwtint_t ref_b,   ref_e;
	int	query_b, query_e;
	int	rid;
//      approximate  value 
	int	n_mis ;
	int	n_gap ;
	int	n_ext ;
// 	approximate  score
	int	app_score ;
	int	sw_score  ;
} aln_res_t ;

typedef struct{
	int m,n ;
	aln_res_t *a ;
} aln_res_v ;

#define aln_chain_t_init(at,size) ((at).m  = (size) , (at).n = 0 , (at).a =  malloc((size)*sizeof(aln_seed_t)))

#define aln_chain_v_init(av)  ( (av).n = 0 , (av).m = 1 , (av).a = malloc(sizeof(aln_chain_t)),  aln_chain_t_init((av).a[0],(av).size) )

#define round_size(x)  ( (x)--, (x) |= (x) >> 1, (x) |= (x) >> 2 , (x) |= (x)>> 4 , (x) |= (x)>> 8 ,  (x)|= (x)>>16 , (x)++)


#define chain_at_add(at,s_aln) do{\
	if((at).m ==  (at).n)  (at).m =  (at).m << 1 ,  (at).a = realloc((at).a , (at).m*sizeof(aln_seed_t)) ; \
	at.a[at.n++] = s_aln ;  \
}while(0)\

#define chain_av_add(av,s_aln) do{\
	(av).n++;\
	if((av).m ==  (av).n)  (av).m =  (av).m << 1 , (av).a = realloc((av).a , (av).m*sizeof(aln_chain_t)) ; \
	aln_chain_t_init((av).a[(av).n], (av).size);\
	chain_at_add((av).a[(av).n],s_aln); \
}while(0)\

#define chain_av_last(av)   (av).a[(av).n].a[(av).a[(av).n].n -1]
#define chain_av_last_c(av,s_aln) (av).a[(av).n].a[(av).a[(av).n].n -1] = s_aln

#define aln_chain_v_free(av)  do{\
				int  chain_i = 0 ;\
				for(;chain_i < (av).n + 1;chain_i++) \
					free((av).a[chain_i].a);\
				free((av).a);\
}while(0)\


int aln_core(const opt_t *opt);
int cm_pe(const opt_t *opt , aln_res_v *rev);
int print_sam();


//  for  unit test  
typedef struct {
	int	len ;
	char	 op ;
}  cigar_t ;

typedef struct{
	int	beg , end ;
	int	right ;
}  key_pos_t ;

typedef struct{
	key_pos_t *val ;
	int    n;
} key_pos_v ;

void  print_av_info(aln_chain_v av , const opt_t *opt);
void  print_at_info(aln_chain_t *at , const opt_t *opt);
int  test_bns(const aln_chain_v av , const seq_t *seq ,const opt_t *opt );
int	test_pos(char *name ,  const aln_chain_v  av , int sel , int l_seed ,const opt_t *opt);


#endif
