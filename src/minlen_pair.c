#include <zlib.h>
#include <stdint.h>
#include <stdio.h>
#include "kseq.h"
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[]){
  gzFile fp, rp, out1, out2;
  kseq_t *seq1, *seq2;
  int l;
  if ( argc != 6 ){
    fprintf(stderr, "Usage: %s <minlen> <in_R1.seq> <in_R2.seq> <out_R1.fq.gz> <out_R2.fq.gz>\n", argv[0]);
    return 1;
  }
  int minlen = atoi(argv[1]);
  fp = gzopen(argv[2], "r"); // STEP 2: open the file handler
  rp = gzopen(argv[3], "r");
  out1 = gzopen(argv[4], "w");
  out2 = gzopen(argv[5], "w");
  seq1 = kseq_init(fp);
  seq2 = kseq_init(rp);
  while ( ((l = kseq_read(seq1)) >= 0) &&
          ((l = kseq_read(seq2)) >= 0)) { // STEP 4: read sequence
    if ( (seq1->seq.l >= minlen) && (seq2->seq.l >= minlen) ){
      //print seq1 out
      gzprintf(out1, "@%s %s\n", seq1->name.s, seq1->comment.s);
      gzprintf(out1, "%s\n", seq1->seq.s);
      gzprintf(out1, "+\n");
      gzprintf(out1, "%s\n", seq1->qual.s);
      //print seq2 out
      gzprintf(out2, "@%s %s\n", seq2->name.s, seq2->comment.s);
      gzprintf(out2, "%s\n", seq2->seq.s);
      gzprintf(out2, "+\n");
      gzprintf(out2, "%s\n", seq2->qual.s);
    }
  }
  kseq_destroy(seq1); // STEP 5: destroy seq
  kseq_destroy(seq2);
  gzclose(fp); // STEP 6: close the file handler
  gzclose(rp);
  gzclose(out1);
  gzclose(out2);
  return 0;
}
