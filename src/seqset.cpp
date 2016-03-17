#include "seqset.h"
#include <math.h>
#include <limits.h>
#include <string.h>


ss_type	readSeqSet(char *filename,a_type A, seq_format format, int MAX_SEQS, const bool addTermini)
{
	Alignment alignment(filename, A, format, addTermini, MAX_SEQS);

	return alignment.getSeqSet();
}

void generateSeqSets(char *filename, a_type A, bool readPvals, seq_format format, int MAX_SEQS, ss_type* posSet, ss_type* negSet){
	fprintf(stderr,  "Option autoThreshold not supported anymore !\n");
	exit(-1);
	//Alignment alignment(filename, A, format);

}

void calculateSequenceFeatures(ss_type set, a_type A){
	set->total[0] = 0;
	double avgLength = 0;
	for(int i=1; i<=set->nent; i++){
		e_type E = set->entity[i];
		int not_X = 0;
		for(int j=1; j<= E->n; j++){
			if(E->S[0][j] != 0) not_X++;
		}
		if(not_X!=0){
		   avgLength += log(not_X);
		}
		set->total[i] = set->total[i-1] + not_X;
		set->avgLength[i] = exp(avgLength/i);
	}
}

void filter_lowComplexity(ss_type set, a_type A){
	int minRepLength = 50;
	int filtered=0;
	int total= 0;
	for(int s=1; s<=set->nent; s++){
		uint8_t* seq = set->entity[s]->S[0];
		int seqlen = set->entity[s]->n;
		int repLength = 2;
		int nucs[2];
		nucs[0] = seq[1];
		nucs[1] = seq[2];
		total += seqlen;
		for(int i=3; i<= seqlen + 1; i++){
			//fprintf(stderr, "%d/%d:%d ", nucs[0], nucs[1], repLength);
			int base;
			if( i== seqlen+1) base = 0;
			else base = seq[i];
			if(nucs[0] == base || nucs[1] == base){
				repLength++;
			}else if(nucs[0] == nucs[1]){
				nucs[1] = base;
				repLength++;
			}else{
				if(repLength > minRepLength && nucs[0] != 0 && nucs[1] != 0){
					//fprintf(stderr, "mask low complexit of length %d in %s (%d, %d)\n",
					//		repLength, set->entity[s]->info[0], nucs[0], nucs[1]);
					for(int j=1; j<=repLength; j++){
						seq[i-j] = 0;
						filtered++;
					}
				}else{
					i -= (repLength-1);
					nucs[0] = seq[i++];
					nucs[1] = seq[i];
				}
				repLength = 2;
			}
		}
//		fprintf(stderr, "\n");
//		for(int i=1; i<=seqlen; i++) fprintf(stderr, "%c", AlphaChar(seq[i], A));
//		fprintf(stderr, "\n");
	}
	std::cerr << filtered << " of " << total << " nucleotides masked because of low complexity" << std::endl;
}

/* filter repeats from the sequence */
void filter_repeats(ss_type set, a_type A)
{
	char* cmask = (char*)malloc((set->max_leng+1)*sizeof(char));
	char* id = (char*)malloc(100*sizeof(char));

	for(int s=1; s<=set->nent; s++){
		uint8_t* seq = set->entity[s]->S[0];
		int seqlen = set->entity[s]->n;

		for(int i=1; i<= seqlen; i++) cmask[i]=0;

		for(int repLength = 20; repLength <= 50; repLength++){

			int W = repLength;
			int MM = W/20;

			if( FilterRep(seq, seqlen, cmask, id, repLength, W, MM) ){
				//fprintf(stderr, "seq: %d, repLength: %d, repeats: %d, W: %d, MM: %d\n", s, repLength, repeats, W, MM);
			}
		}

		// mask sequences
		uint8_t maskedChar = nAlpha(A) == 4 ? AlphaCode('N', A) : AlphaCode('X', A);
		bool masked = false;
		for(int i=1; i<=seqlen; i++)if(cmask[i] != 0) masked = true;

		if(masked){
			//fprintf(stderr, "before masiking: \n");
			//for(int i=1; i<=seqlen; i++) fprintf(stderr, "%c", AlphaChar(seq[i], A));
			//fprintf(stderr, "\n");

			for(int i=1; i<=seqlen; i++){
				if(cmask[i] != 0) seq[i] = maskedChar;
			}

			//fprintf(stderr, "after masiking: \n");
			//for(int i=1; i<=seqlen; i++) fprintf(stderr, "%c", AlphaChar(seq[i], A));
			//fprintf(stderr, "\n\n");
		}
	}
	free(cmask);
	free(id);
}

void createRevcomp(ss_type set, a_type A){
	if(nAlpha(A) != 4){
		fprintf(stderr, "reverse complement only for nucleotides\n");
		exit(-1);
	}
	for(int i=1; i<=set->nent; i++){
		//fprintf(stderr, "%d: ", i);
		int n = set->entity[i]->n;
		set->entity[i]->n = 2*n+1;
		for(int m=0; m<set->entity[i]->mseq; m++){
			/* store memory for reverse complement plus one N position between both strands */
			set->entity[i]->S[m] = (unsigned char*)realloc(set->entity[i]->S[m], n*2 + 2);

			/* set N character in the middle */
			set->entity[i]->S[m][n+1] = 0;

			/* write reverse complement */
			for(int j=n+2, k=n; j<= 2*n+1; j++, k--){
				switch(set->entity[i]->S[m][k]){
					case 0:	set->entity[i]->S[m][j] = 0; break;
					case 1: set->entity[i]->S[m][j] = 4; break;
					case 2: set->entity[i]->S[m][j] = 3; break;
					case 3: set->entity[i]->S[m][j] = 2; break;
					case 4: set->entity[i]->S[m][j] = 1; break;
				}
			}

			//for(int j=1; j<= set->entity[i]->n; j++){
			//	fprintf(stderr, "%c", AlphaChar(set->entity[i]->S[m][j], A));
			//}
			//fprintf(stderr, "\n");
		}
	}
	/* update seq set information */
	set->max_leng = set->max_leng*2 + 1;
	set->min_leng = set->min_leng*2 + 1;

	double avgLength = 0;
	set->total[0] = 0;
	for(int i=1; i<= set->nent; i++){
		e_type seq = set->entity[i];
		int not_X = 0;
		for(int j=1; j<= seq->n; j++){
			if(seq->S[0][j] != '0') not_X++;
		}
		if(not_X!=0) avgLength += log(not_X);

		set->total[i] = set->total[i-1] + not_X;
		set->avgLength[i] = exp(avgLength/i);
	}
}

void fillGapsInMultipleAlignment(ss_type set, a_type A){
	for(int i=1; i<=set->nent; i++){
		e_type ali = set->entity[i];

		for(int m=ali->mseq-2; m>=1; m--){
			for(int pos = 1; pos <= ali->n; pos++){
				/* overwrite all gap positions in not target specie with nucleotides from more distantly related species */
				if(ali->S[m][pos]==0) ali->S[m][pos] = ali->S[m+1][pos];
			}
		}

		for(int m=ali->mseq-1; m>0; m--){
			bool notAligned = true;
			for(int j=1; j<= set->entity[i]->n; j++){
				if(set->entity[i]->S[m][j] != 0){
					notAligned = false;
					break;
				}
			}
			if(notAligned){
				free(ali->S[m]);
				free(ali->info[m]);
				ali->mseq--;
			}
		}
	}
}

void filterMaxMultipleSequences(ss_type set, int max_multiple_sequences){
	set->max_MultSeq = 0;
	for(int i=1; i<=set->nent; i++){
		e_type ali = set->entity[i];
		for(int m = ali->mseq-1; m >= max_multiple_sequences; m--){
			free(ali->S[m]);
			free(ali->info[m]);
			ali->mseq--;
		}
		set->max_MultSeq = std::max((int)ali->mseq, (int)set->max_MultSeq);
	}
}

ss_type	NilSeqSet(ss_type P)
{
	long i;
	for(i=1; i<=P->nent; i++){
		if(P->entity[i] != NULL) NilSeq(P->entity[i]);
	}
	free(P->name);
	free(P->entity);
	free(P->total);
	free(P->avgLength);
	free(P);
	return (ss_type) NULL;
}

e_type	NilSeq(e_type E)
{
	int i=0;
	if(E!=NULL){
		for(i=0;i<E->mseq;i++){
			free(E->info[i]);
			free(E->S[i]);
		}
		free(E->S);
		free(E->info);
		free(E);
	}
	return NULL;
}
