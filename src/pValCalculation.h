#ifndef _PVALCALC_H_
#define _PVALCALC_H_

#include "Globals.h"
#include "LogTable.h"
#include "NullModel.h"
#include "refinementPhase/Motif.h"

class PVal_Calculator{

public:
	static PVal_Calculator &getInstance();

	template<class KGEN> double calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
		int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
		double& pValScore, bool mops);
	template<class KGEN> double __calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
		int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
		double& pValScore, bool mops);

	double calculatePvalCons(int seq, int pos, int firstMotifColumn, const motif_columns_type& sortedColumns);
	double calculateFinalPval(int len, int delta_0, double p);
	static double getCombinedPval_log(double log_prod_pcons, int motifNb);
	double getWeightedPval(double p1, double p2, double w);
	static double getWeightedPval_log(double p1, double p2, double w);

	double getConsCorrection(int length) { return pConsCorrection[length]; }
	double getPosPval(int pos) { return positionalProb[pos]; }

	motif_columns_type getSortedColumns(double** pwm, motif_columns_type& columns);
	void initPositionalProbs(int motifLength, int firstMotifColumn, region r, int offset); /* precalculate positional probabilities */
	//void initPositionalProbs_newLength(int motifLength, int offset);

private:
	PVal_Calculator();
	~PVal_Calculator();
	PVal_Calculator(const PVal_Calculator &);             // intentionally undefined
	PVal_Calculator & operator=(const PVal_Calculator &); // intentionally undefined

	int sse_count_mismatches_gaps(const uint8_t* query_profile, const uint8_t* db_sequence,
												   const int dbseq_length, int query_length, int pos);
	double* positionalProb;
	int** tmpArray;
	int** tmpArray2;
	double* sum_i;
	double* pConsCorrection;
	double* pOverrepCorrection;

	uint8_t* result_ali_free;

public:
	double getOverrepCorrection(const int num) { return pOverrepCorrection[num]; }

};

template<class KGEN>
inline double PVal_Calculator::__calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
		int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
		double& pValScore, bool mops){

	double finalPval, pValPos;
  pValScore = kmerHash.get_kmer_probability_all_positions(kmer);

//  fprintf(stderr, "seq: %d, pos: %d, pValue: %f, size: %d\t", seq, pos, pValScore, length);
//	for(int m=0; m<1; m++){
//	    for(int b= 0; b < length; b++)fprintf(stderr, "%c", AlphaChar(kmer[b], Global::A));
//		fprintf(stderr, "\n");
//	}

	if(pValScore == 1 || pValScore == -1) return pValScore;

	pValScore *= pOverrepCorrection[static_cast<int>(cols.size())];

	if(Global::usePositionalProbs ) {
		pValPos = positionalProb[pos]; 						/* position probability from precalculated values */
		if(mops) {
			/* combine Score and Pos pValue: weight pValPos by 0.5
			 * formula: (p1*p2^w - p1^(1/w)*p2*w) / (1-w) */
			finalPval = (pValScore*sqrt(pValPos) - pValScore*pValScore*pValPos*0.5)*2;
			//double p = pValScore*pValPos;
			//finalPval = p*(1-log(p));
		}else{
			finalPval = calculateFinalPval(len, delta_0, pValScore*pValPos);
		}
		//if(finalPval < 0.05) {
		//	fprintf(stderr, "len: %d, delta_0: %d, pValScore*pValPos: %e\n", len, delta_0, pValScore*pValPos);
		//	fprintf(stderr, "%d/%d, pValPos: %e, pValScore: %e, finalPval: %e\n", seq, pos, pValPos, pValScore, finalPval);
		//}
	} else {
		if(mops){
			finalPval = pValScore;
		}else{
			if (Global::fixedPosition) {
				finalPval = pValScore / NullModel::getProbability("$", 1);
			} else {
				finalPval = 1-pow(1-pValScore,avgLength-length+1);
			}
		}
	}


	return finalPval;
}


template<class KGEN>
inline double PVal_Calculator::calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
			int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
			double& pValScore, bool mops){

    memcpy(kmer, Global::posSet->entity[seq]->S[0] + pos, size);
    return __calculatePval(kmerHash, kmer, cols, seq, pos, len, delta_0, size, length, avgLength,
    			pValScore, mops);
}

inline double PVal_Calculator::calculateFinalPval(int len, int delta_0, double p){
	double plateau = pow(1-p*len/delta_0, delta_0);
	double hill = exp(-p*len*sum_i[delta_0]);
	return 1-plateau*hill;
}

inline double PVal_Calculator::getCombinedPval_log(double log_prod_pcons, int motifNb){

	double result = log_prod_pcons;

	double quotient = 1;
	double sum = 1;

	double overflowConst = 1e100;
	const double log_overflowConst = 230.258509299; // log(overflowConst);
	for(int i=1; i<motifNb; i++){
		quotient *= ( (-log_prod_pcons) / i );
		sum += quotient;

		if(sum > overflowConst){
			sum /= overflowConst;
			quotient /= overflowConst;
			result += log_overflowConst;
		}
		//fprintf(stderr, "i: %d, quotient: %e, sum: %e\n", i, quotient, sum);
	}
	//if(sum <= 0) fprintf(stderr, "log_prod_pcons: %f (%.2e), motifNb: %d, result: %.2e, sum: %.2e, final: %.2e\n", log_prod_pcons, exp(log_prod_pcons), motifNb, result, sum, result+log(sum));
	assert(sum > 0);
	//if(motifNb == 9025)
	//	fprintf(stderr, "result: %e + log(sum): %e = %e (%e)\n", result, log(sum), result+log(sum), exp(result+log(sum)));

	return result + log(sum); // log( p * sum from 0 to n-1 { (-ln p)^i / i! } )
}


inline double PVal_Calculator::getWeightedPval_log(double p1, double p2, double w){
	/* formula: (p1*p2^w - p1^(1/w)*p2*w) / (1-w) */
	double pComb = p1 + w*p2 + log(1.0/(1-w) * (1-w*pow(exp(p1+w*p2),1.0/w-1)));
	return pComb;
}

inline double PVal_Calculator::getWeightedPval(double p1, double p2, double w){
	/* formula: (p1*p2^w - p1^(1/w)*p2*w) / (1-w) */
	double pComb = (p1*pow(p2,w) - pow(p1,1/w)*p2*w) / (1-w);
	fprintf(stderr, "PV: comb(%g, %g) = %g\n", p1, p2, pComb);
	return pComb;
}

inline int PVal_Calculator::sse_count_mismatches_gaps(const uint8_t* query_profile,
											   const uint8_t* db_sequence,
											   const int dbseq_length,
											   int query_length, int pos
											   //unsigned char* results,
											   )
{
	return 0;

}

#endif
