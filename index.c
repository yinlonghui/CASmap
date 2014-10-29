/*
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>

#include "bntseq.h"
#include "utils.h"
#include "bwt.h"
#include "time.h"
/*
 *  First , test read seq
 */
int is_bwt(ubyte_t *T, int n);

int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

bwt_t *bwt_pac2bwt(const char *fn_pac)
{
	bwt_t *bwt;
	ubyte_t *buf, *buf2;
	int i, pac_size;
	FILE *fp;

	// initialization
	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	bwt->seq_len = bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + 15) >> 4;
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
	buf2 = (ubyte_t*)calloc(pac_size, 1);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	memset(bwt->L2, 0, 5 * 4);
	buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
	for (i = 0; i < bwt->seq_len; ++i) {
		buf[i] = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
		++bwt->L2[1+buf[i]];
	}
	for (i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
	free(buf2);

	// Burrows-Wheeler Transform
	bwt->primary = is_bwt(buf, bwt->seq_len);
	bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
	for (i = 0; i < bwt->seq_len; ++i)
		bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
	free(buf);
	return bwt;
}


#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

void bwt_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) {
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) buf[k++] = bwt->bwt[i/16]; // 16 == sizeof(uint32_t)/2
		++c[bwt_B00(bwt, i)];
	}
	// the last element
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + sizeof(bwtint_t) == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;
}

int bwa_bwt2sa(int argc, char *argv[]) // the "bwt2sa" command
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

static int usage(){
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: cm_index <in.fasta>\n");
	fprintf(stderr,"Options: none\n");
	fprintf(stderr,"\n");
	return 0 ;
}
int main( int argc , char *argv[])
{
	if(argc != 2) return usage();
	double  t_start ,t_end ;
	int algo_type;
	char *tmp_name[2];
	init_run_time();
	t_start =  get_run_time();
	
	gzFile fp = xzopen(argv[1],"r");
	tmp_name[0] =  calloc(strlen(argv[1])+15,sizeof(char));
	tmp_name[1] =  calloc(strlen(argv[1])+15,sizeof(char));
	
	/* read  reference to pac */
	fprintf(stderr,"[cm_index] Pack fasta...\n");
	int64_t l_pac =  bns_fasta2bntseq(fp,argv[1],0);
	t_end   =  get_run_time();
	fprintf(stderr,"Finish pack fasta %.2f sec \n",t_end-t_start);
	err_gzclose(fp);
	/* Generate BWT */
	
	t_start =  get_run_time();
	fprintf(stderr, "[cm_index] Construct BWT for the packed sequence...\n");
	strcpy(tmp_name[0],argv[1]);
	strcat(tmp_name[0],".pac");
	strcpy(tmp_name[1],argv[1]);
	strcat(tmp_name[1],".bwt");
	algo_type = l_pac > 50000000 ? 0 : 1;


	if(algo_type){
		bwt_t *bwt;
		bwt = bwt_pac2bwt(tmp_name[0]);
		bwt_dump_bwt(tmp_name[1],bwt);
		bwt_destroy(bwt); 
/*    
 *    Must destory  bwt_bwtgen have difference struct..
 */
	}else  bwt_bwtgen(tmp_name[0],tmp_name[1]); // ..
	t_end   =  get_run_time();
	fprintf(stderr,"Finish Construct BWT %.2f sec \n",t_end-t_start);


	t_start =  get_run_time();
	bwt_t *bwt;
	fprintf(stderr, "[cm_index] Update BWT... ");
	bwt = bwt_restore_bwt(tmp_name[1]);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(tmp_name[1], bwt);
	bwt_destroy(bwt);
	t_end   =  get_run_time();
	fprintf(stderr," %.2f sec \n",t_end-t_start);

	t_start =  get_run_time();
	gzFile fr = xzopen(argv[1], "r");
	fprintf(stderr, "[cm_index] Pack forward-only FASTA... ");
	l_pac = bns_fasta2bntseq(fr, argv[1], 1);
	err_gzclose(fr);
	t_end   =  get_run_time();
	fprintf(stderr," %.2f sec \n",t_end-t_start);

	t_start =  get_run_time();
	strcpy(tmp_name[0], argv[1]); strcat(tmp_name[0], ".sa");
	fprintf(stderr, "[cm_index] Construct SA from BWT and Occ... ");
	bwt = bwt_restore_bwt(tmp_name[1]);
	bwt_cal_sa(bwt, 32);
	bwt_dump_sa(tmp_name[0], bwt);
	bwt_destroy(bwt);
	t_end   =  get_run_time();
	fprintf(stderr," %.2f sec \n",t_end-t_start);

	free(tmp_name[0]);
	free(tmp_name[1]);
	return 0 ;
}
