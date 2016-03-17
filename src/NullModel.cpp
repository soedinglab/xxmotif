#include <cstddef>
#include <iostream>

#include "NullModel.h"
#include "Globals.h"
#include "UngappedKmer.h"

float NullModel::_alpha = 2;
uint8_t NullModel::_asize = 4;
float NullModel::_beta = 0;
int* NullModel::_coeffs = NULL;
float* NullModel::_conds = NULL;
float* NullModel::_counts = NULL;
float NullModel::_factor = 8; /* 4 * 2 as _asize and _alpha default to 4 and 2 resp. */
int NullModel::_fields = 156; /* 5^0 + 5^1 + 5^2 + 5^3 (1 + mono-/di-/trinucleotides) as _order defaults to 2 */
float* NullModel::_freqs = NULL;
bool NullModel::_gaps = true;
int NullModel::_order = 2;
float* NullModel::_probs = NULL;

int NullModel::_bits = 3;
float* NullModel::_rconds = NULL;
int NullModel::_rfields = 512;
float* NullModel::_rprobs = NULL;

void NullModel::expProbs(){

	for( int k=1; k < _fields; k++ ){
		_probs[k] = expf( _probs[k] );
		_conds[k] = expf( _conds[k] );
	}
}
double NullModel::getProbability(char const * const kmer, const int length) {
  static unsigned char trans[100];
  const int len = (length == -1) ? static_cast<int> (strlen(kmer)) : length;
  for (int i = 0; i < len; ++i) {
    trans[i] = AlphaCode((unsigned char)kmer[i], Global::A);
  }
  return getProbability(trans, len);
}

double NullModel::getProbability(unsigned char const * const kmer,
    const int length) {
  double p;
  if (length <= _order + 1) {
    p = _rprobs[UngappedKmer(kmer, length)];
  } else {
    p = _rprobs[UngappedKmer(kmer, _order + 1)];
    for (int stop = _order + 1; stop < length; ++stop) {
      p *= _rconds[UngappedKmer(kmer + stop - _order, _order + 1)];
    }
  }
  return p;
}

void NullModel::logProbs(){

	for( int k=1; k < _fields; k++ ){
		_probs[k] = logf( _probs[k] );
		_conds[k] = logf( _conds[k] );
	}
}

void NullModel::init(ss_type sequences, float alpha, float beta, int order,
    bool gaps, float* freqs) {
  if (alpha < 0) {
    fprintf(stderr, "Negative pseudocounts factor (alpha): %f\n", alpha);
    exit(1);
  }
  _alpha = alpha;

  if (beta < 0) {
    fprintf(stderr, "Negative counts offset (beta): %f\n", beta);
    exit(1);
  }
  _beta = beta;

  if (order < 0) {
    fprintf(stderr, "Negative order: %d\n", order);
    exit(1);
  }
  _order = order;
  _gaps = gaps;
  _freqs = freqs;
  init(sequences);
}

void NullModel::init(ss_type sequences) {
  _asize = static_cast<uint8_t> (nAlpha( Global::A ));
  _factor = static_cast<float> (_asize) * _alpha;
  _bits = static_cast<int> (ceil(log(_asize + 1) / log(2)));
  UngappedKmer::setNumberOfBits(_bits);

  uint8_t* kmer = (uint8_t*) calloc(_order + 1, sizeof(uint8_t));
  uint8_t* new_kmer = (uint8_t*) calloc(_order + 1, sizeof(uint8_t));

  const int base = _asize + _gaps;

  _fields = 1 + base;
  for (int i = 2, powBase = base; i < _order + 2; i++) {
    powBase *= base;
    _fields += powBase;
  }

  _counts = static_cast<float*> (calloc(_fields, sizeof(float)));
  _conds = static_cast<float*> (calloc(_fields, sizeof(float)));
  _probs = static_cast<float*> (calloc(_fields, sizeof(float)));
  _coeffs = static_cast<int*> (calloc(_order + 2, sizeof(int)));

  // _coeffs[k] == pow( base, k );
  for (int k = 0, powBase = 1; k < _order + 2; k++) {
    _coeffs[k] = powBase;
    powBase *= base;
  }

  /* Restructuring parameters */
  _rfields = static_cast<int> (pow(2, (_order + 1) * _bits));
  _rconds = static_cast<float*> (calloc(_rfields, sizeof(float)));
  _rprobs = static_cast<float*> (calloc(_rfields, sizeof(float)));

  /* Calculate counts */
  const int N = sequences->nent;
    for (int n = 1; n <= N; n++) { /* across sequences */
      const int L = sequences->entity[n]->n;
      unsigned char * const s = sequences->entity[n]->S[0];
      for (int l = 1; l <= L; l++) { /* across positions */
        for (int k = 0; k <= _order && l + k <= L && s[l + k]; k++) {
          /* s[l+k]? = char <n>? */
          _counts[sub2ind(s, l, k)]++;
        }
      }
    }

	/* Printouts */
//	printHoNullModel( _counts, _fields-1, offsets );

	/* Substract counts offset ( beta ) */
	if (_beta > 0) {
    for (int i = 1; i < _fields; i++) {
      _counts[i] = std::max(0.0f, (float) (_counts[i] - _beta));
    }
  }

	/* Printouts */
//	printHoNullModel( _counts, _fields-1, offsets );

	/* Calculate total counts for 0th order */
	float ncounts = 0;
  for (int i = 1; i <= _asize; i++) {
    ncounts += _counts[i];
  }

	/* Calculate probs */
	int order = 0; /* 1-mers */
  if (order <= _order) {
    if (_freqs != NULL) {
      for (int n = 1; n <= _asize; n++) {
        _conds[n] = _probs[n] = (_counts[n] + _factor * _freqs[n]) / (ncounts
            + _factor);
      }
    } else {
      _freqs = (float*) calloc(_asize + 1, sizeof(float));
      for (int n = 1; n <= _asize; n++) { /* no pseudocounts in 0th order */
        _conds[n] = _probs[n] = _freqs[n] = _counts[n] / ncounts;
      }
    }
  }

	/* 1-mers++ */
	for (order++; order <= _order; order++) {
    calculateKmerProbs(kmer, 0, order - 1);
  }

	/* Printouts */
  //	printHoNullModel( _conds, _fields-1, offsets );
  //	printHoNullModel( _probs, _fields-1, offsets );

	/* Calculate probs (for gapped kmers) */
	if (_gaps) {
    order = 0; /* 1-mers */
    if (order <= _order) {
      _conds[_asize + 1] = _probs[_asize + 1] = 1;
    }
    /* 1-mers++ */
    for (order++; order <= _order; order++)
      calculateGappedKmerProbs(kmer, 0, order, -1);
  }

	/* Printouts */
//	printHoNullModel( _conds, _fields-1, offsets );
//	printHoNullModel( _probs, _fields-1, offsets );

	/* Copy (conditional) probabilities to Holgers index structure */
	for (order = _order; order >= 0; order--) {
    copyProbs(kmer, new_kmer, 0, order);
  }

	/* Frees */
	free(kmer);
  free(new_kmer);
}

/* Call: calculateGappedKmerProbs( kmer[order+1], 0, order, -1 )*/
void NullModel::calculateGappedKmerProbs( unsigned char* kmer, int pos, int order, int lastGap ){

	if( pos > order ){
		calculateGappedKmerProbs( kmer, order, lastGap );
	}
	else{
		unsigned char k;

		for( k=1; k <= _asize; k++ ){

			kmer[pos] = k;
			calculateGappedKmerProbs( kmer, pos+1, order, lastGap );
		}
		kmer[pos] = k;
		calculateGappedKmerProbs( kmer, pos+1, order, pos );
	}
}

void NullModel::calculateGappedKmerProbs( unsigned char* kmer, int order, int lastGap ){

	int i, ii, k;
	uint8_t X;

	uint8_t N = static_cast<uint8_t>(_asize+1);
	i = sub2ind( kmer, order );

	if( lastGap > -1 )
		_conds[i] = _probs[i] / _probs[sub2ind(kmer,order-1)];

	for( k=lastGap+1; k <= order; k++ ){

		X = kmer[k]; kmer[k] = N;
		ii = sub2ind( kmer, order );

		_probs[ii] += _probs[i];

		kmer[k] = X;
	}
}

/* Call: calculateKmerProbs( kmer[order+1], 0, order-1 )*/
void NullModel::calculateKmerProbs( unsigned char* kmer, int pos, int order_minus_1 ){

	for( uint8_t k=1; k <= _asize; k++ ){

		kmer[pos] = k;

		if( pos == order_minus_1 )
			calculateKmerProbs( kmer, order_minus_1 );
		else
			calculateKmerProbs( kmer, pos+1, order_minus_1 );
	}
}

void NullModel::calculateKmerProbs( unsigned char* kmer, int order_minus_1 ){

	int i, ii, p;
	int order;
	float conds_ii, S;

	int* index = (int*)calloc(_asize+1, sizeof(int));
	int* iindex = (int*)calloc(_asize+1, sizeof(int));

	S = 0;
	order = order_minus_1 + 1;

	for( uint8_t k=1; k <= _asize; k++ ){
		kmer[order] = k;

		ii = sub2ind( kmer, order );
		iindex[k] = ii;

		i = sub2ind( kmer, order-1 );
		index[k] = i;

		p = sub2ind( kmer+1, order-1 );

		conds_ii = ( _counts[ii] + _factor*_conds[p] ) / ( _counts[i] + _factor );
		S += conds_ii;

		_conds[ii] = conds_ii;
	}
	for( uint8_t k=1; k <= _asize; k++ ){

		i = index[k];
		ii = iindex[k];

		_conds[ii] /= S;
		_probs[ii] = _conds[ii] * _probs[i];
	}
	free(index); free(iindex);
}

void NullModel::copyProbs( unsigned char* kmer, unsigned char* new_kmer, int pos, int order ){

	if( pos > order ){
		copyProbs( kmer, new_kmer, order );
	}
	else{
		uint8_t k;

		for( k=1; k <= _asize; k++ ){

			kmer[pos] = k;
			new_kmer[pos] = k;

			copyProbs( kmer, new_kmer, pos+1, order );
		}
		kmer[pos] = k;
		new_kmer[pos] = 0;

		copyProbs( kmer, new_kmer, pos+1, order );
	}
}

void NullModel::copyProbs( unsigned char* kmer, unsigned char* new_kmer, int order ){

	int i = sub2ind( kmer, order );
	int new_i = sub2ind2( new_kmer, order );

	_rprobs[new_i] = _probs[i];
	_rconds[new_i] = _conds[i];
}
