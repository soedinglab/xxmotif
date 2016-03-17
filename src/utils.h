#ifndef UTILS
#define UTILS

#include "Globals.h"
#include "LogTable.h"
#include <cmath>
#include <cstdarg>
#include <string>
#include <execinfo.h>


void getHigherOrderAS(int i, int dim, char* aa);

double log_binomial(const int N, const int k, const double prob, const double* LOG_sum_i);
double binomial(const int N, const int k, const double prob, const double* LOG_sum_i);
double calculate_log_bonferonni(const motif_columns_type& motif_columns, const double log_neff);
double log_binomial_coefficient(const int N, const int k, const double* LOG_sum_i);

double calculateOrderStatisticsPvalue(const int motifNb, const int N, const double p);
double calculateDistanceFromRandom(const int motifNb, const int N, const double p);
double calculateHyperDistribution(int b, int N, int B, int n);

double getPseudocounts(double pseudo, int K);
float fast_log(float x);

/*================================================================================*/


inline void getHigherOrderAS(int i, int dim, char* aa){
	for(int j=dim; j>=1;j--){
		int id = i;
		if(j>1) id /= (int)pow(4, j-1);
		i -= id*(int)pow(4, j-1);
		aa[dim-j] = AlphaChar(id+1,Global::A);
	}
	aa[dim] = '\0';
}

inline double calculate_log_bonferonni(const motif_columns_type& motif_columns, const double log_neff){
		double gapOpening = Global::gapOpening;
		double gapExtension = Global::gapExtension;

		double gapPenalty = 0;

		motif_columns_type::const_iterator it = motif_columns.begin();
		motif_columns_type::const_iterator next = motif_columns.begin(); next++;
		for(; next != motif_columns.end(); it++, next++){
			int gapSize = (*next - *it - 1);
			if(gapSize > 0){
				gapPenalty += gapOpening;
				gapPenalty += (gapSize-1)*gapExtension;
			}
		}
		//fprintf(stderr, "%d * %f = %f (%e)\n", (int)motif_columns.size(), log_neff, motif_columns.size()*log_neff, exp(motif_columns.size()*log_neff));
		return (motif_columns.back() - motif_columns.front() + 1)*log_neff;
}



/* calculate (N over k) * p^k * (1-p)^(N-k) */
inline double log_binomial(const int N, const int k, const double prob, const double* LOG_sum_i){
	//return ( LOG_sum_i[N] - LOG_sum_i[k] - LOG_sum_i[N-k] + k*fast_log(float(prob)) + (N-k)*fast_log(1-(float)prob) );
	return ( LOG_sum_i[N] - LOG_sum_i[k] - LOG_sum_i[N-k] + k*log(prob) + (N-k)*log(1-prob) );
}

/* calculate (N over k) * p^k * (1-p)^(N-k) */
inline double binomial(const int N, const int k, const double prob, const double* LOG_sum_i){
	return exp(log_binomial(N, k, prob, LOG_sum_i));
}

/* calculate log (N over k) */
inline double log_binomial_coefficient(const int N, const int k, const double* LOG_sum_i){
	return LOG_sum_i[N] - LOG_sum_i[k] - LOG_sum_i[N-k];
}

/*************
 * calculates the pvalue of using seqNb of sequences with order statistics
 * t_k+1 = p/(1-p) * (N-k)/(k+1)*t_k <=> tk = t_k+1 ( p/(1-p) * (N-k) / (k+1) )^-1
 *************/
inline double calculateOrderStatisticsPvalue(const int motifNb, const int N, const double p){
	if(N*p > motifNb || p==1) return 0;

//	fprintf(stderr, "motifNb: %d\tN: %d\tp: %g\n", motifNb, N, p);
	double sum_K_1, delta_K_1;
	int K = 0;

	/* calculate pval for finding the motif in the observed number of sequences */
	double logN = log_binomial(N, motifNb, p, LogTable::LOG_sum_i);

	double alpha_K_1,sum_K=0;
	while(K < N - motifNb ){
		alpha_K_1 = (p*(N-(K+motifNb)))/((1-p)*((K+motifNb)+1));
		sum_K_1 = (1 / alpha_K_1 ) * (sum_K + 1);
		sum_K = sum_K_1;
		delta_K_1 = log( 1 + (1/sum_K_1) );
		if(fabs(logN * 1e-4) > delta_K_1 || logN > 0) break;

		logN += delta_K_1;
		K++;
	}
	return logN;
}

/************
 * calculate distence from random
 * **********/
inline double calculateDistanceFromRandom(const int motifNb, const int N, const double p){
	return (motifNb*1.0/N) / p;
}

inline double calculateHyperDistribution(int b, int N, int B, int n){
	double prob = 0;
	for(int i=b; i<=std::min(n, B); i++){
		prob += exp(log_binomial_coefficient(n, i, LogTable::LOG_sum_i) + \
				log_binomial_coefficient(N-n, B-i, LogTable::LOG_sum_i) - \
				log_binomial_coefficient(N, B, LogTable::LOG_sum_i));
	}
	return log(prob);
}

/*
 * This function returns ln with a max abolute deviation of ± 1.5E-5
 * It takes 1.42E-8 s  whereas log2(x) takes 9.5E-7 s. It is hence 9.4 times faster.
 * It makes use of the representation of 4-byte floating point numbers:
 * seee eeee emmm mmmm mmmm mmmm mmmm mmmm
 * s is the sign,
 * the following 8 bits, eee eee e, give the exponent + 127 (in hex: 0x7f).
 * The following 23 bits, m, give the mantisse, the binary digits behind the
 * decimal point.
 * In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
 * The expression (((*(int *)&x) & 0x7f800000 ) >>23 )-0x7f is the exponent eeeeeeee,
 * i.e. the largest integer that is smaller than log2(x) (e.g. -1 for 0.9). *(int *)&x
 * is an integer which
 * contains the bytes as the floating point variable x is represented in memory.
 * Check:  assert( sizeof(f) == sizeof(int) );
 * Check:  assert( sizeof(f) == 4 );
 */
inline float fast_log(float x) {

	if( x != x ) return x; // nan

	static float lg2[1025];   // lg2[i] = log2[1+x/1024]
	static float diff[1025];  // diff[i]= (lg2[i+1]-lg2[i])/8096 (for interpolation)
	static bool initialized;

	if (!initialized) {
		float prev = 0.0f;
		lg2[0] = 0.0f;
		for (int i = 1; i <= 1024; ++i) {
			lg2[i] = (float)log(float(1024+i))*1.442695041f-10.0f;
			diff[i-1] = (lg2[i]-prev)*1.2352E-4f;
			prev = lg2[i];
		}
		initialized=true;
	}

	int32_t a = (((*((int32_t *)&x)) & 0x7F800000) >>23 )-0x7f;
	int32_t b =  ((*((int32_t *)&x)) & 0x007FE000) >>13;
	int32_t c =  ((*((int32_t *)&x)) & 0x00001FFF);

	return ((float)a + lg2[b] + diff[b]*(float)(c)) * 0.69314718f;
}

// sprintf-like helper that returns a string.
inline std::string strprintf(const char* str, ...) {
  char *buffer = new char[1000];
  va_list ap;
  va_start(ap, str);
  vsprintf(buffer, str, ap);
  va_end(ap);
  std::string rv(buffer);
  delete [] buffer;
  return rv;
}

inline std::string translateSeq(const unsigned char* kmer, const int len) {
	std::string s;
	for (int i=0; i<len; ++i) {
		s += AlphaChar(kmer[i], Global::A);
	}
	return s;
}

/* Obtain a backtrace and print it to stdout. */
inline void print_trace (void) {
  void *array[10];
  int size;
  char **strings;
  int i;
  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);
  printf ("Obtained %d stack frames.\n", size);
  for (i = 0; i < size; i++)
     printf ("%s\n", strings[i]);
  free (strings);
}

/* remove all characters in a string from another string */
inline std::string strip_all(const std::string &s, const std::string &chars) {
	std::string ret(s);
	for (std::string::iterator it=ret.begin(); it != ret.end(); ) {
		for (size_t i=0; i<chars.length(); ++i) {
			if (*it == chars[i]) {
				it = ret.erase(it);
				goto CONTINUE_OUTER;
			}
		}
		++it;
		CONTINUE_OUTER: ;
	}
	return ret;
}

inline double getPseudocounts(double ps, int K){
	const double pseudo = ps * K;
	return pseudo;
}

inline std::string strip_spaces(const std::string &s) {
	return strip_all(s, " ");
}

#endif
