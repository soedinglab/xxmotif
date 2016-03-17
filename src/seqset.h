#ifndef SEQSET
#define SEQSET

#include "alphabet.h"
#include "seqFormat/Alignment.h"
#include <stdint.h>

/******************************* PRIVATE *************************************/
void	seqset_error(const char *s);
bool 	FilterRep(uint8_t *seq, int seqlen, char* cmask, char* id, int p, int W, int MM);
e_type	NilSeq(e_type E);			/* undefine entity */

/******************************* PUBLIC *************************************/
ss_type  readSeqSet(char *filename,a_type A, seq_format format, int MAX_SEQS, const bool addTermini = false);
void 	 generateSeqSets(char *filename,a_type A, bool readPvals, seq_format format, int MAX_SEQS, ss_type* posSet, ss_type* negSet);

void 	 calculateSequenceFeatures(ss_type set, a_type A);
void 	 filter_repeats(ss_type set, a_type A);
void     filter_lowComplexity(ss_type set, a_type A);
void 	 createRevcomp(ss_type set, a_type A);
void 	 fillGapsInMultipleAlignment(ss_type set, a_type A);
void 	 filterMaxMultipleSequences(ss_type set, int max_multiple_sequences);
long     *LengthsSeqSet(ss_type P);

ss_type  NilSeqSet(ss_type P);


/* Filter short repeats of periodicity p with max MM mismatches out of W
   pairs compared (i.e. window length = W+p) */
inline bool FilterRep(uint8_t *seq, int seqlen, char* cmask, char* id, int p, int W, int MM)
{
  bool masked = false;
  int i, j, k;
  int sum;         // number of identities

  for (j=0; j<W; j++) id[j]=0;
  j=sum=0;

  for (i=1+p; i<=seqlen; i++){
      sum -= id[j];
      id[j] = ( seq[i-p]==seq[i]? 1 : 0 );
      sum += id[j];
      if (sum >= W-MM){
    	  for (k=i-W-p+1; k<=i; k++){
    		  masked = true;
    		  if(k>0) cmask[k]=1;
    	  }
      }
      if (++j>=W){
    	  j=0;
      }
  }
  return masked;
}

#define	LenSeq(E)		((E)->n)
#endif
