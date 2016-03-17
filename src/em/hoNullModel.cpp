#include "hoNullModel.h"

float   hoNullModel::_alpha = 10.0f;
int*	hoNullModel::_coeffs = NULL;
double*	hoNullModel::_conds	= NULL;
double*	hoNullModel::_counts = NULL;
double*	hoNullModel::_countsx = NULL;
int		hoNullModel::_fields = 5; // 1 + nAlpha( Global::A )
double*	hoNullModel::_freqs = NULL;
int*	hoNullModel::_offsets = NULL;
int		hoNullModel::_order = 0;
double*	hoNullModel::_probs = NULL;

void hoNullModel::expProbs(){

	for( int k=1; k < _fields; k++ ){
		_probs[k] = exp( _probs[k] );
		_conds[k] = exp( _conds[k] );
	}
}

void hoNullModel::logProbs(){

	for( int k=1; k < _fields; k++ ){
		_probs[k] = log( _probs[k] );
		_conds[k] = log( _conds[k] );
	}
}

void hoNullModel::init( ss_type sequences, float alpha, int order, double* freqs ){

	if( alpha < 0 ){
		fprintf( stderr, "Background pseudocounts factor < 0" );
		exit(0);
	}
	_alpha = alpha;

	if( order < 0 ){
		fprintf( stderr, "Background order < 0" );
		exit(0);
	}
	_order = order;

	_freqs = freqs;

	init( sequences );
}

void hoNullModel::init( ss_type sequences ){

	int i, k;
	int l, L;
	int n, N;

	int order;
	double ncounts;

	unsigned char* s;

	uint8_t *kmer = new uint8_t[_order+1];

	_coeffs = ( int* )calloc( _order+2, sizeof( int ) );
	_offsets = ( int* )calloc( _order+2, sizeof( int ) );

	for( k=0; k < _order+2; k++ ){
		_coeffs[k] = static_cast<int>( pow( nAlpha( Global::A ), k ) );
		_offsets[k] = k ? _offsets[k-1]+_coeffs[k] : _coeffs[k];
	}
	_fields = _offsets[k-1];

	_counts = ( double* )calloc( _fields, sizeof( double ) );
	_countsx = ( double* )calloc( _offsets[_order], sizeof( double ) );

	_conds = ( double* )calloc( _fields, sizeof( double ) );
	_probs = ( double* )calloc( _fields, sizeof( double ) );

	/*
	 * calculate counts
	 */
	N = sequences->nent;
	for( n=1; n <= N; n++ ){ // across sequences

		s = sequences->entity[n]->S[0];
		L = sequences->entity[n]->n;

		if( _order >= L ){
			fprintf( stderr, "Background sequence no. %d shorter than null "
					 "model order %d\n", n, _order );
			exit(-1);
		}

		for( l=1; l <= L; l++ ){ // across positions
			for( k=0; k <= _order && l+k <= L && s[l+k]; k++ ){
				// sequence[l+k]
				// checks whether nucleotide equals N
				// represented by 0
				_counts[sub2ind( s, l, k )]++;
				if( k > 0 ){
					// kmer not followed by N character
					_countsx[sub2ind( s, l, k-1 )]++;
				}
			}
			if( k <= _order && l+k > L ){
				// last column handling
				_countsx[sub2ind( s, l, k-1 )]++;
			}
		}
	}

	/*
	 * calculate 0th-order probabilities
	 */
	order = 0;
	if( order <= _order ){

		ncounts = 0.0;
		for( i=1; i <= nAlpha( Global::A ); i++ ){
			ncounts += _counts[i];
		}

		if( _freqs != NULL ){
			for( n=1; n <= nAlpha( Global::A ); n++ ){
				_conds[n] = _probs[n] = ( _counts[n] + _alpha*_freqs[n] ) /
						                ( ncounts + _alpha );
			}
		}
		else{
			_freqs = ( double* )calloc( nAlpha( Global::A )+1, sizeof( double ) );
			for( n=1; n <= nAlpha( Global::A ); n++ ){ // no pseudocounts
				_conds[n] = _probs[n] = _freqs[n] = _counts[n] / ncounts;
			}
		}

		order++;
	}

	/*
	 * calculate higher-order probabilities
	 */
	for( ; order <= _order; order++ ){
		calculateKmerProbs( kmer, 0, order );
	}

	if( Global::verbose ){
		print();
	}

	delete[] kmer;
}

void hoNullModel::save( char* baseFileName ){

	/**
	 * Save interpolated Markov background model to three flat files
	 * (1) baseFileName.freqsBg (background frequencies)
	 * (2) baseFileName.probsBg (probabilities)
	 * (3) baseFileName.condsBg (conditional probabilities)
	 */

	std::stringstream str;

	str.str( "" );
	str << Global::outputDirectory << "/" << baseFileName << ".freqsBg";
	FILE* f_frequencies = fopen( str.str().c_str(), "w" );

	str.str( "" );
	str << Global::outputDirectory << "/" << baseFileName << ".condsBg";
	FILE* f_conditionals = fopen( str.str().c_str(), "w" );

	str.str( "" );
	str << Global::outputDirectory << "/" << baseFileName << ".probsBg";
	FILE* f_probabilities = fopen( str.str().c_str(), "w" );

	int i, k;

	for( i=1; i <= nAlpha( Global::A ); ++i ){
		fprintf( f_frequencies, "%f ", _freqs[i] );
	}
	fprintf( f_frequencies, "\n" );
	fclose( f_frequencies );

	for( k=1, i=1; i < _fields; ++i ){
		if( i==_offsets[k] ){
			fprintf( f_conditionals, "\n" );
			fprintf( f_probabilities, "\n" );
			++k;
		}
		fprintf( f_conditionals, "%f ", _conds[i] );
		fprintf( f_probabilities, "%f ", _probs[i] );
	}
	fprintf( f_conditionals, "\n\n" );
	fprintf( f_probabilities, "\n\n" );

	fclose( f_conditionals );
	fclose( f_probabilities );
}

void hoNullModel::calculateKmerProbs( unsigned char* kmer, int pos, int order ){

	for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
		kmer[pos] = k; // next kmer base
		if( pos == ( order - 1 ) ){
			// calculate probabilities
			calculateKmerProbs( kmer, order );
		} else{
			// add another base
			calculateKmerProbs( kmer, pos+1, order );
		}
	}
}

void hoNullModel::calculateKmerProbs( unsigned char* kmer, int order ){

	int i, ii, p;

	int *indices_i = new int[nAlpha( Global::A )+1];
	int *indices_ii = new int[nAlpha( Global::A )+1];

	double normFactor = 0.0;
	for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
		// last kmer base
		kmer[order] = k;
		// shorter kmer indices
		i = sub2ind( kmer, order-1 );
		// longer kmer indices
		ii = sub2ind( kmer, order );
		// pseudo-counts kmer indices
		p = sub2ind( kmer+1, order-1 );

		indices_ii[k] = ii;
		indices_i[k] = i;

		// calculate conditional probabilities
		_conds[ii] = ( _counts[ii] + _alpha*_conds[p] ) /
				     ( _countsx[i] + _alpha );
		normFactor += _conds[ii];
	}
	for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){

		i = indices_i[k];
		ii = indices_ii[k];

		// normalize conditional probabilities
		_conds[ii] /= normFactor;
		// calculate probabilities
		_probs[ii] = _conds[ii] * _probs[i];
	}

	delete[] indices_i;
	delete[] indices_ii;
}

void hoNullModel::print(){
	printf( " ____________\n"
			"|            |\n"
			"| NULL MODEL |\n"
			"|____________|\n" );
	printf( " ________\n"
			"|        |\n"
			"| COUNTS |\n"
			"|________|\n\n" );
	printInterpolatedMarkovBackgroundModel( _counts, _fields-1, _offsets );
	printf( " _________\n"
			"|         |\n"
			"| COUNTSX |\n"
			"|_________|\n\n" );
	printInterpolatedMarkovBackgroundModel( _countsx, _offsets[_order]-1,
			                                _offsets );
	printf( " ___________________________\n"
			"|                           |\n"
			"| CONDITIONAL PROBABILITIES |\n"
			"|___________________________|\n\n" );
	printInterpolatedMarkovBackgroundModel( _conds, _fields-1, _offsets );
	printf( " _____________________\n"
			"|                     |\n"
			"| TOTAL PROBABILITIES |\n"
			"|_____________________|\n\n" );
	printInterpolatedMarkovBackgroundModel( _probs, _fields-1, _offsets );
	printf( "\n" );
}
