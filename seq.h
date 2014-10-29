#ifndef __SEQ_H
#define __SEQ_H



#include <stdio.h>
#include <stdint.h>

typedef uint32_t  seq_sz_t;
typedef struct {
	char	*seq;
	char	*name;
	char	*add;
	char	*qual;
	size_t	alen;
	size_t	qlen;
	size_t	nlen;
	size_t  len;
} seq_t;

typedef void *seqdb_t;

seqdb_t seq_open(const char *seq_fn);
int seqt_release(seq_t *seq);
void seqdb_get(seqdb_t db,seq_t **pseq);
int seqdb_close(seqdb_t db);
int seq_reverse(int len ,unsigned char *seq);
int seq_complement(int len ,unsigned char *seq);
int char2nt4(unsigned char **byte,seq_sz_t  len , const char *seq,int set_al);

#define xeorror(msg) fprintf(stderr,"%s:%s:%s:%d\n",msg,__FUNCTION__,__FILE__,__LINE__)

#endif
