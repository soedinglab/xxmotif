#ifndef NULLMODEL_H_
#define NULLMODEL_H_

#include <stdlib.h>
#include "seqFormat/Alignment.h"

class NullModel{

public:

	static double getProbability(unsigned char const * const kmer, const int length);
	static double getProbability(char const * const kmer, const int length = -1);
	static double getConditional(uint64_t index){ return _rconds[index]; }
	static int getOrder() { return _order; }

	static void destruct(){
		free( _coeffs ); free( _conds ); free( _counts ); free( _freqs ); free( _probs ); free( _rconds ); free( _rprobs );
	}

	static void expProbs();

	static float* getConds(){
		return _conds;
	}

	static float* getCounts(){
		return _counts;
	}

	static float* getFreqs(){
		return _freqs;
	}

	static float* getProbs(){
		return _probs;
	}

	static float getAlpha() {
	  return _alpha;
	}

	static float getBeta() {
	  return _beta;
	}

	static void init( ss_type sequences, float alpha, float beta, int order, bool gaps, float* freqs );

	static void init( ss_type sequences );

	static void logProbs();

	/* Calculates array indices for k-mers */
	static int sub2ind( unsigned char const * const sequence, int order, bool byrow=true );

	/* Calculates array indices for sequence k-mers */
	static int sub2ind( unsigned char const * const sequence, int pos, int order, bool byrow=true );

private:

	/* Calculates (conditional) probabilities for gapped kmers */
	static void calculateGappedKmerProbs( unsigned char* kmer, int pos, int order, int lastGap );
	static void calculateGappedKmerProbs( unsigned char* kmer, int order, int lastGap );

	/* Calculates (conditional) probabilities for kmers */
	static void calculateKmerProbs( unsigned char* kmer, int pos, int order_minus_1 );
	static void calculateKmerProbs( unsigned char* kmer, int order_minus_1 );

	/* Restructure (conditional) probabilities structure */
	static void copyProbs( unsigned char* kmer, unsigned char* new_kmer, int pos, int order );
	static void copyProbs( unsigned char* kmer, unsigned char* new_kmer, int order );

	/* Calculates Holger's array indices for k-mers */
	static int sub2ind2( unsigned char* sequence, int order );

	/* Pseudocounts factor */
	static float _alpha;

	/* Alphabet size */
	static uint8_t _asize;

	/* Counts offset */
	static float _beta;

	/* Address coefficients */
	static int* _coeffs;

	/* HO kmer conditionals */
	static float* _conds;

	/* HO kmer counts */
	static float* _counts;

	/* Pseuodcounts factor */
	static float _factor;

	/* Field number in _conds/_counts/_probs arrays */
	static int _fields;

	/* Monomer background frequencies */
	static float* _freqs;

	/* Calculate statistics for kmers with gaps? */
	static bool _gaps;

	/* Model order */
	static int _order;

	/* HO kmer probabilities */
	static float* _probs;

	static int _bits;
	static float* _rconds;
	static int _rfields;
	static float* _rprobs;
};

inline int NullModel::sub2ind( unsigned char const * const sequence, int pos, int order, bool byrow ){

	int i, k, l;
	i = 0;

	if( byrow ){
		for( k=0, l=order; k<=order; ++k, --l ){
			i += _coeffs[l] * sequence[pos+k];
		}
	}
	else{
		for( k=0; k<=order; ++k ){
			i += _coeffs[k] * sequence[pos+k];
		}
	}

	return i;
}

inline int NullModel::sub2ind( unsigned char const * const sequence, int order, bool byrow ){

	int i, k, l;
	i = 0;

	if( byrow ){
		for( k=0, l=order; k<=order; ++k, --l ){
			i += _coeffs[l] * sequence[k];
		}
	}
	else{
		for( k=0; k<=order; ++k ){
			i += _coeffs[k] * sequence[k];
		}
	}

	return i;
}

inline int NullModel::sub2ind2( unsigned char* sequence, int order ){

	int res = *sequence; /* kmer[0] */

	for( int i=1; i <= order; ++i ){

		res <<= _bits;
		res += *(sequence + i); /* kmer[i] */
	}

	return res;
}

#endif /* HONULLMODEL_H_ */
