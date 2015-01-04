#ifndef __ALN_H
#define __ALN_H
#include "seq.h"
#include "bwa.h"
#include "kvec.h"
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
	/* record seqdb and  refdb */
	file_seq_t *fs;
	file_ref_t *fr;
	/*  match(a) and mismatche(b) score  , and matrix */
	int8_t mat[25];
	int  a , b;
	/*  penalty gap and extension*/
	int  o_del , o_ins ;
	int  e_ins , e_del ;
	int  overlap ;
	/*  extension local alignment */
	int  pen_clip5 , pen_clip3 ;
	int  zdrop ;
	int  w;
	/*   nunmber of need reads  */
	int  n_need ;
	/* length seed (k-mer)*/
	int  l_seed ;
	int  verbose;
	int  flag;
} opt_t ;


/*
 *  fastq structure 
 */

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
	int	len ;
} aln_seed_t;


#define  is_extend(at) ((at).flags & ALLOW_EXTEND)
#define  is_overlap(at) ((at).flags & HAVE_OVERLAP)
#define  DISALLOW_EXTEND   (~ALLOW_EXTEND)

typedef struct {
	aln_seed_t *a  ;
	int n , m ;
} aln_seed_v ;

typedef struct{
	bwtint_t ref_b,   ref_e;
	int	query_b, query_e;
	int	rid;
//      approximate  information 
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

typedef struct {
	aln_seed_t  *a ;
	int	flags;
	int	seedcov;
	int	m,n;
} aln_chain_t ;

#define  ALLOW_EXTEND  1 
#define  HAVE_OVERLAP  2

typedef struct {
	aln_chain_t *a;
	int     m,n,size,offset;
} aln_chain_v;

//   The marco handle chain vector

#define aln_chain_t_init(at,size) ((at).m  = (size) , (at).n = 0 , (at).a =  malloc((size)*sizeof(aln_seed_t)))

#define aln_chain_v_init(av)  ( (av).n = 0 , (av).m = 1 , (av).a = malloc(sizeof(aln_chain_t)),  aln_chain_t_init((av).a[0],(av).size) )

#define round_size(x)  ( (x)--, (x) |= (x) >> 1, (x) |= (x) >> 2 , (x) |= (x)>> 4 , (x) |= (x)>> 8 ,  (x)|= (x)>>16 , (x)++)
//  too  slow ..
#define get_top_value(v)  ((v).a[(v).n-1])

#define push_chain_seed(av,seed,type,i)  kv_push(type,(av).a[i],seed)

#define push_chain(av,v_tpye, seed,s_tpye) do{\
	(av).n++;\
	if((av).n == (av).m ){\
		(av).m =  (av).m << 1 ;  \
		(av).a = realloc((av).a ,sizeof(v_tpye)*(av).m);\
	}\
	aln_chain_t_init((av).a[(av).n],(av).size); \
	kv_push(s_tpye,(av).a[(av).n],seed);\
}while(0)



#define aln_chain_v_free(av)  do{\
				int  chain_i = 0 ;\
				for(;chain_i < (av).n + 1;chain_i++) \
					free((av).a[chain_i].a);\
				free((av).a);\
}while(0)\





int aln_core(const opt_t *opt);
int cm_pe(const opt_t *opt , aln_res_v *rev[2]);
int print_sam();


//  import  other file function 



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

//   for unit test fuction 
extern void  print_av_info(aln_chain_v av , const opt_t *opt);
extern void  print_at_info(aln_chain_t *at);
extern void  print_res_info(aln_res_t  *at);
extern void  print_resv_info(aln_res_v *av);
extern int  test_bns(const aln_chain_v av , const seq_t *seq ,const opt_t *opt );
extern int	test_pos(char *name ,  const aln_chain_v  av , int sel , int l_seed ,const opt_t *opt);
extern void  unit_sv_seed(aln_seed_v *kv_seed ,const opt_t *opt);
extern int	unit_extend_1( char *name , aln_res_v  *rev , int sel , int l_seed , const opt_t *opt);
extern int	find_mismatch( char *name , int l_seed , int sel);
extern int	find_insertion( char *name , int l_seed , int sel);
extern int	find_deletion( char *name , int l_seed , int sel);

#endif
