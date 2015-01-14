#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "seq.h"


char nt4_table[32]={
	0,4,1,4,4,4,2,4,
	4,4,4,4,4,4,4,4,
	4,4,4,3,4,4,4,4,
	4,4,4,4,4,4,4,4,
};
#define a_to_nt4(c) (((c)>'a')?nt4_table[(c)-'a']:nt4_table[(c)-'A'])

typedef	struct {
	FILE  *fin ;
}seqdb;

int char2nt4(unsigned char **byte, seq_sz_t len , const char *seq, int set_al)
{
	if(len <= 0) 
	return -1;
	seq_sz_t i ;
	*byte = malloc(sizeof(unsigned char)*len);
	unsigned char *tmp;
	for( i = 0 , tmp = *byte ; i < len ; i++ , tmp = tmp + 1){
		*tmp = a_to_nt4(*(seq+i));
	}
	return 0 ;
}
#define a_to_a_r(seq)  (seq > 3 ? seq: 3 -seq)
int seq_reverse(int len ,unsigned char *seq)
{
	int i = 0;
	for( i = 0 ; i < (len>>1) + len%2  ; i++){
		char tmp = a_to_a_r(seq[len-1-i]); 
		seq[len-1-i] = a_to_a_r(seq[i]);
		seq[i] = tmp;

	}

	return 0;
}
int seq_complement(int len ,unsigned char *seq)
{
	int i ;
	for( i = 0 ; i < len ; i++)
		seq[i] = a_to_a_r(seq[i]);
	

	return 0;
}
int seqdb_close(seqdb_t db)
{
	seqdb *pdb = db ;
	if( pdb->fin == NULL){
		return 1 ;
	}
	fclose(pdb->fin);
	free(pdb);
	return 0 ;
}
static seqdb *_seq_open(const char *seq_fn)
{
	seqdb *pdb = NULL;
	pdb = malloc(sizeof(seqdb));
	pdb->fin = fopen(seq_fn,"r");
	if(pdb->fin == NULL){
		xeorror("open");
		goto ERR_0;
	}
	return pdb;
ERR_0:
	free(pdb);
	return NULL;
}
seqdb_t seq_open(const char *seq_fn)
{
	seqdb *pdb = NULL;
	if(seq_fn){
		pdb = _seq_open(seq_fn);
		if(pdb == NULL){
			seqdb_close(pdb);
			return NULL;
		}

	}
	return pdb;
}

void _seqdb_get(seqdb *db,seq_t **p)

{
	seq_t *pseq ;
	size_t len ;
	ssize_t rc;
	*p = calloc(1,sizeof(seq_t));
	pseq = *p ;
	if(pseq == NULL)
		xeorror("malloc");
	if(db->fin == NULL){
		xeorror("fin is NULL");
		goto ERR0;
	}
	rc = getline(&pseq->name,&pseq->nlen,db->fin);
	pseq->nlen = rc-1;
	if(rc <= 0)
		goto ERR0;
	rc = getline(&pseq->seq,&pseq->len,db->fin);
	pseq->len = rc-1;
	if(rc <= 0)
		goto ERR0;
	rc = getline(&pseq->add,&len,db->fin);
	pseq->alen = rc-1;
	if(rc<=0)
		goto ERR0;
	rc = getline(&pseq->qual,&len,db->fin);
	pseq->qlen = rc-1;
	if(rc<=0)
		goto ERR0;
	return ;
ERR0:
	free(pseq->name);
	free(pseq->seq);
	free(pseq->qual);
	free(pseq->add);
	free(pseq);
	*p = NULL;
}

void  seqdb_get(seqdb_t db ,seq_t **pseq)
{
	seqdb *pdb = db;
	_seqdb_get(pdb,pseq);
}
int seqt_release(seq_t *seq)
{
	if(seq){
		free(seq->seq);
		free(seq->name);
		free(seq->sam);
		free(seq->add);
		free(seq->qual);
		free(seq);
	}
	return 0;
}
