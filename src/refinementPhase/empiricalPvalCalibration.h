#ifndef EMPIRICAL_PVAL_CALIBRATION
#define EMPIRICAL_PVAL_CALIBRATION

#include "Motif.h"
#include "StartPosUpdater.h"
#include "../Globals.h"
#include <limits>

class empiricalPvalCalibration{
public :
	static empiricalPvalCalibration& getInstance() {
		static empiricalPvalCalibration instance;
		return instance;
	}

	void reinitialize(Motif* motif, unsigned char* kmer);
	double get_kmer_probability_all_positions(unsigned char* kmer);

private:
	empiricalPvalCalibration();
	~empiricalPvalCalibration();
	empiricalPvalCalibration(const empiricalPvalCalibration &);             // intentionally undefined
	empiricalPvalCalibration & operator=(const empiricalPvalCalibration &); // intentionally undefined

	void updateRefScoreList(float smax);
	double checkMotifOnRefScoreList(float smax);
	float getMatrixScore(unsigned char* kmer);

	Motif* M;
	sorted_sites* scoreList;				/* scores used for empirical calibration */
	double* bgLog;
	float smax;
	int refSetNb;
};

inline void empiricalPvalCalibration::reinitialize(Motif* motif, unsigned char* kmer){
	M = motif;
	bgLog = Global::negSet == NULL ? Global::posBg_log : Global::negBg_log;
	smax = std::numeric_limits<float>::max();
	smax = getMatrixScore(kmer) + 1e-5f; /* avoid precision errors in float */
	updateRefScoreList(smax);
}

inline double empiricalPvalCalibration::get_kmer_probability_all_positions(unsigned char* kmer){
	return checkMotifOnRefScoreList(getMatrixScore(kmer));
}

inline float empiricalPvalCalibration::getMatrixScore(unsigned char* kmer){
	double score = 0;
	for(int pos=0; pos < M->getMotifLength(); pos++){
	    score -= M->getPWM()[pos+M->getFirstMotifColumn()][kmer[pos]] - bgLog[kmer[pos]];
	}
	if(score > smax){
		fprintf(stderr, "\nERROR: %f > %f\n", score, smax);
		exit(-1);
	}
	return static_cast<float>(score);
}

#endif
