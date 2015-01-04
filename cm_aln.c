#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>

#include "aln.h"
#include "seq.h"
#include "time.h"


static int usage(){

	fprintf(stderr,"\n");
	fprintf(stderr,"usgae: cm_aln  <in.fasta>  <in.fastq> [in2.fastq]\n");
	fprintf(stderr,"Options: -f FILE  file to wirte output to instead of stdout \n");
//	fprintf(stderr,"         -n  INT  number of mismatch (seeding) \n");
	fprintf(stderr,"         -k  INT  length of seed(k-mer)\n");

	fprintf(stderr,"         -A  INT  score for a sequence match\n");
	fprintf(stderr,"         -B  INT  penalty for a mismatch \n");
	fprintf(stderr,"         -I  INT  insertion gap open penalty \n");
	fprintf(stderr,"         -i  INT  insertion gap extend penalty \n");
	fprintf(stderr,"         -D  INT  deletion gap open penalty \n"); 
	fprintf(stderr,"         -d  INT  deletion gap extend penalty \n");
	fprintf(stderr,"         -w  INT   band width for banded alignment \n");

	fprintf(stderr,"         -v  INT  output debug info. \n ");
	fprintf(stderr,"                  1、 test only one read for SE alignment.\n");
	fprintf(stderr,"                      print all stack info to trace.\n");
	fprintf(stderr,"\n");
	return  EXIT_FAILURE ;
}


static void  gen_mat(int a , int b, int8_t mat[25])
{
	int i , j  ;
	for( i = 0 ; i < 5 ; i++)
		for( j = 0 ; j < 5 ; j++)
			mat[i*5+j] = (i == 4 ||j == 4) ? -1 :(i!=j ?-b :a);

}
/*
 *   Initialize  Options Structure
 */
static inline opt_t *init_opt()
{
	opt_t *opt =  calloc(1,sizeof(opt_t));
	opt->n_need =  0x40000;
	opt->l_seed = 22 ; //  this seed parameter  need to experiment repeatedly
	opt->a = 1 ;
	opt->b = 4 ;
	opt->pen_clip5 = opt->pen_clip3 =  5 ;
	opt->zdrop =  100 ;
	opt->o_del = opt->o_ins = 6 ;
	opt->e_del = opt->e_ins = 1 ;
	gen_mat(opt->a,opt->b,opt->mat);
	opt->w = 100;
	return opt;
}

int	main(int argc , char *argv[])
{
	opt_t  *opt = init_opt();
/*
 *    parse the command.. 
 */
	int c =  0 , rc = 0 ;
	while( ( c = getopt(argc , argv ,"f:v:z:k:A:B:I:i:D:d:w:")) != -1){
		switch(c){
			/* All state... */
			case 'f':  freopen(optarg,"w",stdout);  break;   
			case 'v':  opt->verbose = atoi(optarg); break;
			case 'z':  opt->n_need = atoi(optarg);  break;
			case 'k':  opt->l_seed = atoi(optarg);  break;

			case 'A':  opt->a =  atoi(optarg); break;
			case 'B':  opt->b =  atoi(optarg); break;
			case 'i':  opt->e_ins =  atoi(optarg); break;
			case 'I':  opt->o_ins =  atoi(optarg); break;
			case 'D':  opt->o_del =  atoi(optarg); break;
			case 'd':  opt->e_del =  atoi(optarg); break;
			case 'w':  opt->w =  atoi(optarg); break;
		
		}
	}
//  test  input  paramter 
#if  0
	fprintf(stderr,"match score  =%d\n", opt->a);
	fprintf(stderr,"mismatch score  =%d\n", opt->b);
	fprintf(stderr,"band  alignment =%d\n",opt->w);
	fprintf(stderr,"gap open   insert  penalty =%d\n", opt->o_ins);
	fprintf(stderr,"gap open   delete  penalty =%d\n", opt->o_del);
	fprintf(stderr,"gap extend insert  penalty =%d\n", opt->e_ins);
	fprintf(stderr,"gap extend insert  penalty =%d\n", opt->e_del);

#endif

// score scheme 
#if 0
	int  i , j , k = 0 ;
	for( i = 0 ; i < 5 ; i++){
		for( j = 0 ; j < 5 ; j++)
			fprintf(stderr,"%d\t",opt->mat[k++]);
		fprintf(stderr,"\n");
	}

#endif
	if(optind + 3 == argc) opt->flag = F_PE  ;
	

	if(optind + 4 > argc && argc > optind + 1){
		double t_start,t_end;
		init_run_time();
		
		t_start = get_run_time();

		//  open  prefix sequence ...
		fprintf(stderr,"Loading bwt and pac.. ");

		opt->fr = malloc(sizeof(file_ref_t));
		opt->fr->fn =  argv[optind];
		opt->fr->idx = bwa_idx_load(argv[optind],BWA_IDX_ALL);
		t_end = get_run_time();

		fprintf(stderr,"%.2f sec.\n",t_end-t_start);

		t_start =  get_run_time();
		
		opt->fs = calloc(opt->flag&F_PE?2:1,sizeof(file_seq_t));
		//  open  read 1 sequence ...
		opt->fs->fn =  argv[optind+1] ; 
		opt->fs->seqdb = seq_open(argv[optind+1]);
	
		if(opt->flag& F_PE){
			// open  read 2  sequence ...
			opt->fs[1].fn = argv[optind+2];
			opt->fs[1].seqdb =  seq_open(argv[optind+2]);
		}
		assert(2*opt->fr->idx->bns->l_pac ==  opt->fr->idx->bwt->seq_len);
		rc = aln_core(opt);
		if(rc == EXIT_FAILURE) goto ALN_FAILURE_EXIT ;
		t_end = get_run_time();
		fprintf(stderr,"Finish the alignment.. total time is  %2.f sec. \n",t_end-t_start);

	}else  goto FAILURE_EXIT; 

//      free opt_t 
	seqdb_close(opt->fs->seqdb);
	if(opt->flag&F_PE) seqdb_close(opt->fs[1].seqdb);
	free(opt->fs);
	bwa_idx_destroy(opt->fr->idx);
	free(opt->fr);
	free(opt);
	return  EXIT_SUCCESS;

FAILURE_EXIT:		
	free(opt);
	return  usage();

ALN_FAILURE_EXIT:

	seqdb_close(opt->fs->seqdb);
	if(opt->flag&F_PE) seqdb_close(opt->fs[1].seqdb);
	free(opt->fs);
	
	bwa_idx_destroy(opt->fr->idx);
	free(opt->fr);
	free(opt);
	return EXIT_FAILURE;
}

