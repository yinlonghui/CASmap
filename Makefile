CC = gcc
DEBUG = -g  
CFLAGS?= $(DEBUG) -Wall -D_FILE_OFFSET_BITS=64 -pthread  -pg
LIBS =  -lm -lz -lpthread -lc
BIN =  cm_index cm_aln 
index_obj = index.o  utils.o bntseq.o bwt_gen.o QSufSort.o bwt.o is.o time.o  
aln_obj = cm_aln.o  aln.o seq.o bwa.o ksw.o kstring.o  utils.o bntseq.o bwt.o malloc_wrap.o time.o
all:$(BIN)

cm_index:$(index_obj)
	$(CC) -o  $@ $^ $(LIBS)	-pg

cm_aln:$(aln_obj)
	$(CC) -o  $@ $^ $(LIBS)	-pg
clean:
	rm $(BIN) *o  *~
