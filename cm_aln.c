#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>

#include "aln.h"
#include "seq.h"
#include "time.h"


static int usage(opt_t *opt){

	fprintf(stderr,"\n");
	fprintf(stderr,"usgae: cm_aln  <in.fasta>  <in.fastq> [in2.fastq]\n");
	fprintf(stderr,"Options: -f FILE  file to wirte output to instead of stdout \n");
	fprintf(stderr,"         -k  INT  length of seed(k-mer)			[default:%d]\n",opt->l_seed);

	fprintf(stderr,"         -A  INT  score for a sequence match		[default:%d]\n",opt->a);
	fprintf(stderr,"         -B  INT  penalty for a mismatch		[default:%d]\n",opt->b);
	fprintf(stderr,"         -I  INT  insertion gap open penalty		[default:%d]\n",opt->o_ins);
	fprintf(stderr,"         -i  INT  insertion gap extend penalty		[default:%d]\n",opt->e_ins);
	fprintf(stderr,"         -D  INT  deletion gap open penalty		[default:%d]\n",opt->o_del); 
	fprintf(stderr,"         -d  INT  deletion gap extend penalty		[default:%d]\n",opt->e_del);
	fprintf(stderr,"	 -e  INT  minimun extend length			[default:%d]\n",opt->min_extend_len);

	fprintf(stderr,"         -w  INT  band width for banded alignment	[default:%d]\n", opt->w);
	fprintf(stderr,"	 -T  INT  threshold score for outputing sam	[default:%d]\n" ,opt->threshold);
	fprintf(stderr,"	 -t  INT  number of threads\n");

	fprintf(stderr,"         -v  INT  output debug info. \n ");
	fprintf(stderr,"\n");
	return  EXIT_FAILURE ;
}


static void  gen_mat(int a , int b, int8_t mat[25])
{
	int i , j  ;
	for( i = 0 ; i < 5 ; i++)
		for( j = 0 ; j < 5 ; j++)
			mat[i*5+j] = (i == 4 ||j == 4) ? -b :(i!=j ?-b :a);

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
	opt->min_extend_len = 30 ;
	gen_mat(opt->a,opt->b,opt->mat);
	opt->w = 100;
	opt->threshold = 30 ;
	opt->n_thread =  1 ;
	return opt;
}

void  cm_print_sam_hdr(const bntseq_t *bns)
{
	int  i ;
	for( i = 0 ; i < bns->n_seqs ; i++){
		printf("@SQ\tSN:%s\tLN:%d\n",bns->anns[i].name,bns->anns[i].len);
	}
}

double t_start,t_end;

int	main(int argc , char *argv[])
{
	opt_t  *opt = init_opt();
/*
 *    parse the command.. 
 */
	int c =  0 , rc = 0 ;
	while( ( c = getopt(argc , argv ,"f:v:z:k:A:B:I:T:i:D:d:w:t:e:")) != -1){
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
			case 'e':  opt->min_extend_len = atoi(optarg); break ;

			case 'T':  opt->threshold =  atoi(optarg) ; break ;
			case 't':  opt->n_thread =  atoi(optarg) ; break ;

			default:   goto  FAILURE_EXIT;	
				   
		
		}
	}

	if(optind + 3 == argc) opt->flag = F_PE  ;
	

	if(optind + 4 > argc && argc > optind + 1){
		init_run_time();
		
		t_start = get_run_time();

		//  open  prefix sequence ...
		fprintf(stderr,"Loading bwt and pac.. ");

		opt->fr = malloc(sizeof(file_ref_t));
		opt->fr->fn =  argv[optind];
		opt->fr->idx = bwa_idx_load(argv[optind],BWA_IDX_ALL);
		t_end = get_run_time();

		fprintf(stderr,"%.2f sec.\n",t_end-t_start);

		cm_print_sam_hdr(opt->fr->idx->bns);
		printf("@PG\tID:CASmap\tPN:CASmap\tVN:0.1\tCL:%s",argv[0]);
		int  i  ;
		for( i  = 1 ; i < argc ; i++ )
			printf(" %s",argv[i]);
		printf("\n");

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
		fprintf(stderr,"Finish the alignment.. total time is  %.2f sec. \n",t_end-t_start);

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
	 
	usage(opt);
	free(opt);
	return EXIT_FAILURE ;

ALN_FAILURE_EXIT:

	seqdb_close(opt->fs->seqdb);
	if(opt->flag&F_PE) seqdb_close(opt->fs[1].seqdb);
	free(opt->fs);
	
	bwa_idx_destroy(opt->fr->idx);
	free(opt->fr);
	free(opt);
	return EXIT_FAILURE;
}

