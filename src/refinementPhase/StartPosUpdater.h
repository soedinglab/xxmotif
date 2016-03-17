#ifndef STARTPOS_UPDATER_H
#define STARTPOS_UPDATER_H

#include "Motif.h"
#include "StartPosUtils.h"
#include "Sorted_Sites.h"
#include "../Globals.h"
#include "../pValCalculation.h"
#include "../ThresholdChecker.h"
#include "../prodpval_stat.h"
#include "../refinementPhase/Motif.h"

#include "../branch_and_bound.h"

class StartPosUpdater
{
public:
  static StartPosUpdater &getInstance();


  void getProbabilitiesPerSite(list<double> &siteProb, list<double> &siteCons, double& combined_pValCons, Motif* m);
  void fill_startPosList_with_bestSubset(Motif* motif, double& bestPval, int& bestMotifNb);
  double fill_startPosList_with_elongSites(StartPosContainer& realStartPosList, int bestMotifNb, motif_columns_type& tmpColumns);
  void set_BindingSites(Motif* motif, const ThresholdChecker &pValThreshold);

  void update_PWM_with_bestSubset(Motif* motif, int& bestMotifNb, double& bestPval, bool plus, double pseudo);
  void update_elongPWM_with_bestSubset(Motif* motif, double **elongPWM, motif_columns_type& tmpColumns,
		  int& bestMotifNb, double& bestPval, bool plus, double pseudo);

  void setPositionalPval(Motif* motif);

  template<class KGEN>
  void fill_sites_with_startPosList(KGEN& kmerHash, Motif* motif, bool instanceThreshold=false);

  template<class KGEN>
  void fill_sites_with_possible_startPositions(KGEN& kmerHash, bool keepPosSetSize, Motif* motif);

  template<class KGEN>
  void fill_elongSites_with_possible_startPositions(KGEN& kmerHash, Motif* motif, double** pwm, motif_columns_type& tmpColumns, bool palin);

  template<class KGEN>
  void fill_sites_with_possible_startPositions_oops(KGEN& kmerHash, Motif* motif);

  void update_sites_with_empirical_Pvals(Motif* motif);

  sorted_sites* getBestStartPos(){ return bestStartPos; }
  sorted_sites* getSortedSites(){ return sortedSites; }
  int* getMotifsCountPerSequence(){ return motifsPerSequenceCount; }


private:
  StartPosUpdater();
  ~StartPosUpdater();
  StartPosUpdater(const StartPosUpdater &);             // intentionally undefined
  StartPosUpdater & operator=(const StartPosUpdater &); // intentionally undefined

  void getBestSubset(double **pwm, motif_columns_type& columns, sorted_sites* sites, int nbSites, int posSetSize,
		  double& bestPval, int& bestMotifNb, int &bestCounter, int& bestMotifNbPlus, int& bestCounterPlus, int oldMotifNb);
  int getBestSubset_plus(double **pwm, motif_columns_type& columns, sorted_sites* sites,
		  int sitesNb, int posSetSize, int oldMotifNb);

  void fill_startPosList_with_sites(StartPosContainer& sites, int firstMotifColumn, int nb);
  void fill_startPosList_with_matches(MatchContainer& sites, int firstMotifColumn, int nb);
  void updateElongPWM(double** pwm, double** elongPWM, motif_columns_type& tmpColumns);

  static int sort(const void *x, const void *y) {
    return (*(double*)x - *(double*)y > 0);
  }

  PVal_Calculator& pValCalculator;

  sorted_sites* sortedSites;        /* stores sites sorted by combined pValues */
  sorted_sites *sortedSites_elong;  /* stores temporal sites used in motif elongation */
  sorted_sites* bestStartPos;		/* stores sites sorted by pValues without positional information */
  double*	pvalCons;				/* stores conservation pValues */
  int* motifsPerSequenceCount;		/* stores how many motifs are found per sequence */
  unsigned char* kmer;

  int sitesAlloc1;		/* have to be stored to free memory at the end of program when seqset already freed */
  int sitesAlloc2;

  int nbSites;
  int nbBestSites;
  int maxSites;
  int refSetNb;
  double log_posSetSize;

  double* log_posSetSites;
};

inline void StartPosUpdater::getProbabilitiesPerSite(list<double> &siteProb, list<double> &siteCons, double& combined_pValCons, Motif* m){
	motif_columns_type sortedColumns = pValCalculator.getSortedColumns(m->getPWM(), m->getMotifColumns());
	double pConsCorrection = pow(Global::consCorrection, m->getMotifLength());

	int region = m->getEnrichment().endRegion - m->getEnrichment().startRegion + 1;
	int length = m->getMotifLength();
	int size = static_cast<int>(sizeof(unsigned char)*length);

	/* init positional probabilities */
	if(Global::usePositionalProbs){
		pValCalculator.initPositionalProbs(m->getMotifLength(), m->getFirstMotifColumn(), m->getEnrichment(), 0);
	}

	double log_prod_pcons = 0;

	Kmer_Generator_Reloaded<aa_hash_t> &kmerHash = Kmer_Generator_Reloaded<aa_hash_t>::getInstance();
	kmerHash.reinitialize(m->getPWM(), m->getMotifColumns(), m->getPosSetSize(), Global::posSet->total[m->getPosSetSize()]);

	for (StartPosContainer::const_iterator it_startPos= m->getStartPosList().begin();
		it_startPos != m->getStartPosList().end(); it_startPos++){

		int pos = it_startPos->pos + m->getFirstMotifColumn() - 1;

		double pValCons = pValCalculator.calculatePvalCons(it_startPos->seq, pos, m->getFirstMotifColumn(), sortedColumns);
		pValCons *= pConsCorrection;
		if(pValCons > 1) pValCons = 1;

		log_prod_pcons += log(pValCons);

		siteCons.push_back(pValCons);

		int len = Global::posSet->entity[it_startPos->seq]->n - length + 1;
		int delta_0 = std::min(len, region);

		double pValScore;
		double finalPval = pValCalculator.calculatePval(kmerHash, kmer, m->getMotifColumns(), \
						it_startPos->seq, pos, len, delta_0, size, length, Global::posSet->avgLength[m->getPosSetSize()], \
						pValScore, Global::multipleOccurrence);
		siteProb.push_back(finalPval);
	}
	combined_pValCons = pValCalculator.getCombinedPval_log(log_prod_pcons, (int)m->getStartPosList().size());
}


inline void StartPosUpdater::getBestSubset(double **pwm, motif_columns_type& columns, sorted_sites* sites, int nbSites, int posSetSize,
		double& bestPval, int& bestMotifNb, int& bestCounter, int& bestMotifNbPlus, int& bestCounterPlus, int oldMotifNb){
	bool debug = false;
	bool plot = false;

	double maxChange = oldMotifNb * 0.3;
	if(maxChange == 0)maxChange = maxSites;

	if(nbSites == 0){
		bestMotifNb = 0;
		return;
	}
	bestMotifNb = 1;
	bestCounter = 1;
	motif_columns_type sortedColumns = pValCalculator.getSortedColumns(pwm, columns);
	int firstMotifColumn = columns.front();
	double log_prod_pcons = 0;
	double log_ranks_pValue = 0;

	double pConsCorrection = pow(Global::consCorrection, static_cast<int>(columns.size()));

	/* reset motifs per sequence counter */
	memset(motifsPerSequenceCount, 0, (Global::posSet->nent+1)* sizeof(int));
	int motifCounter = 0;

	bestPval = 1;

	for(int i=0; i<std::min(nbSites, posSetSize*Global::maxSeqCount); i++){
		if(sites[i]->score == 1) break;

		motifsPerSequenceCount[sites[i]->sequence]++;
		if(motifsPerSequenceCount[sites[i]->sequence] > Global::maxMotifsPerSequence) continue;

		motifCounter++;
		//fprintf(stderr, "%d %d/%d/%f\t", motifCounter, sites[i]->sequence, sites[i]->startPos, sites[i]->score);

		if(Global::posSet->max_MultSeq > 1){
			double pValCons = pValCalculator.calculatePvalCons(sites[i]->sequence, sites[i]->startPos, firstMotifColumn, sortedColumns);
			pValCons *= pConsCorrection;
			if(pValCons > 1) pValCons = 1;
			//if(sites[i]->sequence == 6)fprintf(stderr, "seq: %d, sp: %d, %f\n", sites[i]->sequence, sites[i]->startPos, pValCons);
			//pvalCons[motifCounter-1] = pValCons;
			log_prod_pcons += log(pValCons);
		}

		if(Global::useRankPvalues){
			log_ranks_pValue += log_posSetSites[sites[i]->sequence];
		}
		/* if next site has the same score, this will automatically have a better pValue => skip current site */
		if(i+1< nbSites && sites[i+1]->score == sites[i]->score) continue;

		double pValue = calculateOrderStatisticsPvalue(motifCounter, maxSites, sites[i]->score);

		if(debug){
			for(unsigned int j=0; j< columns.size(); j++) cerr << AlphaChar(Global::posSet->entity[sites[i]->sequence]->S[0][sites[i]->startPos+j], Global::A);
			fprintf(stderr, "\tnb: %d:, maxSites: %d, score: %e => pValue: %f\n", motifCounter, maxSites, sites[i]->score, pValue);
		}
		if(plot){
			fprintf(stderr, "%d\t%.2e\t%.2f\t%.2f\n", motifCounter, sites[i]->score, pValue, 0.0);
		}
		if(Global::useRankPvalues){
			double pValRanks = pValCalculator.getCombinedPval_log(log_ranks_pValue-(i+1)*log_posSetSize, i+1);
			pValue += 0.5*pValRanks + LogTable::LOG_i[2];
			if(debug)fprintf(stderr, "\tpValRanks: %f (%e)", pValRanks, exp(pValRanks));
		}

		if(Global::posSet->max_MultSeq > 1){
			double pValCons;
			//if(Global::useAliFree) pValCons = getCombinedPval(log_prod_pcons, 2*(i+1));
			pValCons = pValCalculator.getCombinedPval_log(log_prod_pcons, i+1);
			/* combined pValue */
			if(debug)fprintf(stderr, "\toverRepPval: %f (%e), sites: %d", pValue, exp(pValue), motifCounter);
			pValue = pValCalculator.getWeightedPval_log(pValue, pValCons, Global::consPvalWeight);
			if(debug)fprintf(stderr, "\tpValCons: %f (%e) \n\tpValCombined = %f (%e)\n\n", pValCons, exp(pValCons), pValue, exp(pValue));
		}


		double plusFrac = 1 + Global::plusFrac;
		if(pValue < bestPval || Global::oneOccurrence){
			bestPval = pValue;
			bestMotifNbPlus = bestMotifNb = motifCounter;
			bestCounterPlus = bestCounter = i+1;
		}else if(plusFrac != 1){
			if(motifCounter > plusFrac*bestMotifNb){
				bestMotifNbPlus = std::max(bestMotifNb+1, static_cast<int>(bestMotifNb*plusFrac + 0.5));
				bestCounterPlus = std::max(bestCounter+1, static_cast<int>(bestCounter*plusFrac + 0.5));
			}else{
				bestMotifNbPlus = motifCounter;
				bestCounterPlus = i+1;
			}
		}
	}



	if(debug){
		fprintf(stderr, "\nbestMotifNb: %d, bestPval: %f (%e), bestMotifNbPlus: %d\n", bestMotifNb, bestPval, exp(bestPval), bestMotifNbPlus);
		exit(-1);
	}
	if(plot){
		fprintf(stderr, "%d\t%.2f\t%.2f\t%.2f\n",maxSites, 1.0, 0.0, 0.0);
		exit(-1);
	}
}


inline void StartPosUpdater::fill_startPosList_with_sites(StartPosContainer& sites, int firstMotifColumn, int nb){
	sites.clear();

	/* reset motifs per sequence counter */
	memset(motifsPerSequenceCount, 0, (Global::posSet->nent+1)* sizeof(int));

	for(int i=0, motifCounter = 0; motifCounter < nb;i++){
		motifsPerSequenceCount[sortedSites[i]->sequence]++;
		if(motifsPerSequenceCount[sortedSites[i]->sequence] > Global::maxMotifsPerSequence) continue;

		++motifCounter;

		sites.push_back(StartPos(sortedSites[i]->sequence, static_cast<int32_t>(sortedSites[i]->startPos-firstMotifColumn+1)));
	}
}

inline void StartPosUpdater::fill_startPosList_with_matches(MatchContainer& sites, int firstMotifColumn, int nb){
	sites.clear();

	for(int i=0; i < nb; i++){
		sites.push_back(Match(sortedSites[i]->sequence, static_cast<int32_t>(sortedSites[i]->startPos-firstMotifColumn+1), sortedSites[i]->score));
	}
}

inline void StartPosUpdater::updateElongPWM(double** pwm, double** elongPWM, motif_columns_type& tmpColumns){
	/* copy current pwm to elong pwm */
	for(motif_columns_type::iterator it_col = tmpColumns.begin(); it_col != tmpColumns.end(); it_col++){
		//cerr << *it_col << " ";
		for(int j = 1; j<= nAlpha(Global::A); j++){
			elongPWM[*it_col][j] = pwm[*it_col][j];
		}
	}
}
#endif
