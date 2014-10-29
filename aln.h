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

int aln_core(opt_t *opt);


#endif
