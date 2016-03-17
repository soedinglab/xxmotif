#ifndef HONULLMODEL_H
#define HONULLMODEL_H

#include "../Globals.h"
#include "../seqFormat/Alignment.h"

#include "hoUtils.h"

class hoNullModel{

public:

	static void destruct(){
		free( _coeffs );
		free( _conds );
		free( _counts );
		free( _countsx );
		free( _freqs );
		free( _offsets );
		free( _probs );
	}

	static void expProbs();

	static int* getCoeffs(){
		return _coeffs;
	}

	static double* getConds(){
		return _conds;
	}

	static double* getFreqs(){
		return _freqs;
	}

	static int* getOffsets(){
		return _offsets;
	}

	static int getOrder(){
		return _order;
	}

	static double* getProbs(){
		return _probs;
	}

	static void init( ss_type sequences, float alpha, int order, double* freqs );

	static void init( ss_type sequences );

	static void logProbs();

	static void save( char* baseFileName );

	/*
	 * calculates oligomer indices
	 */
	static int sub2ind( unsigned char* sequence, int pos, int order, bool byrow=true );

private:

	/*
	 * calculates higher-order probabilities for oligomers
	 */
	static void calculateKmerProbs( unsigned char* kmer, int pos, int order_minus_1 );
	static void calculateKmerProbs( unsigned char* kmer, int order_minus_1 );

	static void print();

	/*
	 * calculates oligomers indices
	 */
	static int sub2ind( unsigned char* sequence, int order, bool byrow=true );

	/*
	 * pseudocounts factor
	 */
	static float _alpha;

	/*
	 * oligomer address coefficients
	 */
	static int* _coeffs;

	/*
	 * oligomer conditional probabilities
	 */
	static double* _conds;

	/*
	 * oligomer counts
	 */
	static double* _counts;

	/*
	 * oligomer counts
	 * counts oligomers without succeeding N character and gaps
	 */
	static double* _countsx;

	/*
	 * oligomer number
	 */
	static int _fields;

	/*
	 * background frequencies
	 * not in use
	 * defaults to frequencies
	 */
	static double* _freqs;

	/*
	 * oligomer address offsets
	 */
	static int* _offsets;

	/*
	 * model order
	 */
	static int _order;

	/*
	 * oligomer probabilities
	 */
	static double* _probs;
};

inline int hoNullModel::sub2ind( unsigned char* sequence, int pos, int order, bool byrow ){

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

inline int hoNullModel::sub2ind( unsigned char* sequence, int order, bool byrow ){

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

#endif /* HONULLMODEL_H */
