CC = gcc
DEBUG = -g
CFLAGS?= $(DEBUG) -Wall -D_FILE_OFFSET_BITS=64 -pthread  
LIBS =  -lm -lz -lpthread 
BIN =  cm_index cm_aln 
index_obj = index.o  utils.o bntseq.o bwt_gen.o QSufSort.o bwt.o is.o time.o  
aln_obj = cm_aln.o  aln.o seq.o bwa.o ksw.o kstring.o  utils.o bntseq.o bwt.o malloc_wrap.o time.o alnpe.o  test.o kthread.o
all:$(BIN)

cm_index:$(index_obj)
	$(CC) -o  $@ $^ $(LIBS)	

cm_aln:$(aln_obj)
	$(CC) -o  $@ $^ $(LIBS)	
clean:
	rm $(BIN) *o  *~
