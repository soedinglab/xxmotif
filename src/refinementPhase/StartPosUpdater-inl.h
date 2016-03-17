#ifndef STARTPOS_UPDATER_INL_H
#define STARTPOS_UPDATER_INL_H

#include "StartPosUpdater.h"

template<class KGEN>
void StartPosUpdater::fill_sites_with_startPosList(KGEN& kmerHash, Motif* motif, bool instanceThreshold){
	double finalPval, pValScore;

	/* init variables */
	int region = motif->getEnrichment().endRegion - motif->getEnrichment().startRegion + 1;
	int length = motif->getMotifLength();
	int size = static_cast<int>(sizeof(unsigned char)*length);

	/* init positional probabilities */
	if(Global::usePositionalProbs){
		pValCalculator.initPositionalProbs(motif->getMotifLength(), motif->getFirstMotifColumn(), motif->getEnrichment(), 0);
	}

	/* if oops model reset probs */
	if(!Global::multipleOccurrence){
		for(int i=0; i<motif->getPosSetSize(); i++){
			sortedSites[i]->score = 1;
		}
	}
	//bool debug = false;
	//if(motif->getTotalSites() == 18 && motif->getPval() < -7.6 && motif->getPval() > -7.7) debug = true;
	//if(debug) cerr << *motif;

	int sites_counter=0;
	for (StartPosContainer::const_iterator it = motif->getStartPosList().begin(); it!= motif->getStartPosList().end(); it++){
		const int32_t seq = it->seq;
		const int32_t startPos = static_cast<int32_t>(it->pos + motif->getFirstMotifColumn() - 1);

		/* initialize sortedSites */
		const int len = Global::posSet->entity[seq]->n - length + 1;
		const int delta_0 = std::min(len, region);

		finalPval = pValCalculator.calculatePval(kmerHash, kmer, motif->getMotifColumns(), \
				seq, startPos, len, delta_0, size, length, Global::posSet->avgLength[motif->getPosSetSize()], \
				pValScore, Global::multipleOccurrence || instanceThreshold);

		if(finalPval >= 1) continue;

		if(Global::multipleOccurrence){
			if(!updateMultOcc(sortedSites, sites_counter, seq, startPos, length, finalPval, Global::posSet)) continue;
			sites_counter++;
		}else if(finalPval < sortedSites[seq-1]->score){
			if(sortedSites[seq-1]->score == 1) sites_counter++;

			sortedSites[seq-1]->sequence = seq;
			sortedSites[seq-1]->startPos = startPos;
			sortedSites[seq-1]->score = static_cast<float>(finalPval);
		}
	}

	if(Global::multipleOccurrence){
		nbSites = sites_counter;
		maxSites = Global::posSet->total[motif->getPosSetSize()] - motif->getPosSetSize()*(length-1);
		if(Global::revcomp){
			maxSites = Global::posSet->total[motif->getPosSetSize()] - 2*motif->getPosSetSize()*(length-1);
		}
	}else{
		nbSites = motif->getPosSetSize();
		maxSites = motif->getPosSetSize();
	}
	log_posSetSize = LogTable::LOG_i[Global::posSet->total[motif->getPosSetSize()]];
}

template<class KGEN>
void StartPosUpdater::fill_sites_with_possible_startPositions(KGEN& kmerHash, bool keepPosSetSize, Motif* motif){
	double finalPval, pValScore;

	/* init variables */
	int region = motif->getEnrichment().endRegion - motif->getEnrichment().startRegion + 1;
	int length = motif->getMotifLength();
	int size = static_cast<int>(sizeof(unsigned char)*length);
	//motif_columns_type sortedColumns = pValCalculator.getSortedColumns(motif->getPWM(), motif->getMotifColumns());

	/* init positional probabilities */
	if(Global::usePositionalProbs){
		pValCalculator.initPositionalProbs(motif->getMotifLength(), motif->getFirstMotifColumn(), motif->getEnrichment(), 0);
	}

	int sites_counter = 0;
	int bestSites_counter = 0;
	int32_t bestPosSetSize = motif->getPosSetSize();

	int32_t sequences = keepPosSetSize ? bestPosSetSize : Global::posSet->nent;
	for(int32_t i=1; i <= sequences; i++){
		/* initialize sites */
		int len = Global::posSet->entity[i]->n - length + 1;
		int delta_0 = std::min(len, region);

		/* test every possible start position */
		for(int32_t j=1; j <= len; j++){
			finalPval = pValCalculator.calculatePval(kmerHash, kmer, motif->getMotifColumns(), \
					i, j, len, delta_0, size, length, Global::posSet->avgLength[motif->getPosSetSize()], \
					pValScore, true);

			if(finalPval >= 1) continue;

			bestStartPos[bestSites_counter]->sequence = i;
			bestStartPos[bestSites_counter]->startPos = j;
			bestStartPos[bestSites_counter]->score = static_cast<float>(pValScore);
			bestSites_counter++;

			if(!updateMultOcc(sortedSites, sites_counter, i, j, length, finalPval, Global::posSet)){
				continue;
			}
			sites_counter++;
		}
	}

	bestStartPos[bestSites_counter]->score = 1;
	while(bestSites_counter > 0 && bestStartPos[bestSites_counter-1]->sequence > bestPosSetSize){ bestStartPos[--bestSites_counter]->score = 1; }
	sortedSites[sites_counter]->score = 1;
	while(sites_counter > 0 && sortedSites[sites_counter-1]->sequence > bestPosSetSize){ sortedSites[--sites_counter]->score = 1; }

	//	for(int i=0; i<sites_counter; i++){
	//		fprintf(stderr, "%d\t", i+1);
	//		for(int j=0; j<motif->getMotifLength(); j++){
	//			fprintf(stderr, "%c", AlphaChar(Global::posSet->entity[sortedSites[i]->sequence]->S[0][sortedSites[i]->startPos+j], Global::A));
	//		}
	//		fprintf(stderr, "\t%d/%d \t %f\n", sortedSites[i]->sequence, sortedSites[i]->startPos, sortedSites[i]->score);
	//	}

	nbSites = std::max(0,sites_counter);
	nbBestSites = std::max(0, bestSites_counter);

	sortSitesByScore(bestStartPos, nbBestSites);

	motif->setPosSetSize(bestPosSetSize);

	maxSites = Global::posSet->total[motif->getPosSetSize()] - motif->getPosSetSize()*(length-1);
	if(Global::revcomp){
		maxSites = Global::posSet->total[motif->getPosSetSize()] - 2*motif->getPosSetSize()*(length-1);
	}
	log_posSetSize = LogTable::LOG_i[Global::posSet->total[motif->getPosSetSize()]];
}


template<class KGEN>
void StartPosUpdater::fill_elongSites_with_possible_startPositions(KGEN& kmerHash, Motif* motif, double** pwm, motif_columns_type& tmpColumns, bool palin){
	double finalPval, pValScore;

	/* init variables */
	int region = motif->getEnrichment().endRegion - motif->getEnrichment().startRegion + 1;
	int length = tmpColumns.back()-tmpColumns.front()+1;
	int size = static_cast<int>(sizeof(unsigned char)*length);
	int firstMotifColumn = motif->getFirstMotifColumn();

	/* init positional probabilities */
	if(Global::usePositionalProbs){
		pValCalculator.initPositionalProbs(tmpColumns.back()-tmpColumns.front()+1, firstMotifColumn, motif->getEnrichment(), tmpColumns.front() - firstMotifColumn);
	}

	/* if oops model reset probs */
	if(!Global::multipleOccurrence){
		for(int32_t i=0; i<motif->getPosSetSize(); i++) sortedSites_elong[i]->score = 1;
	}

	int sites_counter = 0;

	sortSitesByStartPos(bestStartPos, nbBestSites, palin);

	for(int i=0; i<nbBestSites; i++){
		/* initialize sites */
		const int32_t seq = bestStartPos[i]->sequence;
		int32_t pos = static_cast<int32_t>(bestStartPos[i]->startPos + tmpColumns.front() - firstMotifColumn);
		if(palin){
			pos = static_cast<int32_t>(Global::posSet->entity[seq]->n - pos + 1 - length + 1);
		}

		//fprintf(stderr, "%d/%d: ", seq, pos);
		//for(int i=pos; i<pos+length; i++){
		//	fprintf(stderr, "%c", AlphaChar(Global::posSet->entity[seq]->S[0][i], Global::A));
		//}
		//fprintf(stderr, "\n");

		if(pos < 1 || pos + tmpColumns.back() - tmpColumns.front() > Global::posSet->entity[seq]->n) continue;

		int len = Global::posSet->entity[seq]->n - length + 1;
		int delta_0 = std::min(len, region);

		finalPval = pValCalculator.calculatePval(kmerHash, kmer, tmpColumns, \
				seq, pos, len, delta_0, size, length, Global::posSet->avgLength[motif->getPosSetSize()], \
				pValScore, Global::multipleOccurrence);

		if(finalPval >= 1) continue;

		//fprintf(stderr, "final Pval: %f\n", finalPval);

		if(Global::multipleOccurrence){
			if(!updateMultOcc(sortedSites_elong, sites_counter, seq, pos, length, finalPval, Global::posSet)){
				continue;
			}
		}else if(finalPval < sortedSites_elong[seq-1]->score){
			//cerr << i << "\t" << j << "\t" << pValScore << "\t" << Global::posSet->avgLength-length+1 << "\t" << finalPval << endl;
			sortedSites_elong[seq-1]->sequence = seq;
			sortedSites_elong[seq-1]->startPos = pos;
			sortedSites_elong[seq-1]->score = static_cast<float>(finalPval);
		}
		sites_counter++;
	}
	if(Global::multipleOccurrence){
		sortedSites_elong[sites_counter]->score = 1;

		nbSites = sites_counter;
		maxSites = Global::posSet->total[motif->getPosSetSize()] - motif->getPosSetSize()*(length-1);

		if(Global::revcomp){
			maxSites = Global::posSet->total[motif->getPosSetSize()] - 2*motif->getPosSetSize()*(length-1);
		}
	}else{
		nbSites = motif->getPosSetSize();
		maxSites = motif->getPosSetSize();
	}

	log_posSetSize = LogTable::LOG_i[Global::posSet->total[motif->getPosSetSize()]];
}


template<class KGEN>
void StartPosUpdater::fill_sites_with_possible_startPositions_oops(KGEN& kmerHash, Motif* motif){
	double finalPval, pValScore;

	/* init variables */
	int region = motif->getEnrichment().endRegion - motif->getEnrichment().startRegion + 1;
	int length = motif->getMotifLength();
	int size = static_cast<int>(sizeof(unsigned char)*length);
	//motif_columns_type sortedColumns = pValCalculator.getSortedColumns(motif->getPWM(), motif->getMotifColumns());

	/* init positional probabilities */
	if(Global::usePositionalProbs){
		pValCalculator.initPositionalProbs(motif->getMotifLength(), motif->getFirstMotifColumn(), motif->getEnrichment(), 0);
	}
	for(int i=1; i <= Global::posSet->nent; i++){
		sortedSites[i-1]->score = 1;
		sortedSites[i-1]->sequence=i;
		sortedSites[i-1]->startPos=1;

		/* initialize sites */
		const int len = Global::posSet->entity[i]->n - length + 1;
		const int delta_0 = std::min(len, region);

		/* test every possible start position */
		for(int32_t j=1; j <= len; j++){
			finalPval = pValCalculator.calculatePval(kmerHash, kmer, motif->getMotifColumns(), \
					i, j, len, delta_0, size, length, Global::posSet->avgLength[motif->getPosSetSize()], \
					pValScore, false);

			if(finalPval < sortedSites[i-1]->score){
				sortedSites[i-1]->sequence = i;
				sortedSites[i-1]->startPos = j;
				sortedSites[i-1]->score = static_cast<float>(finalPval);
			}
		}
	}

	nbSites = motif->getPosSetSize();
	maxSites = Global::posSet->total[motif->getPosSetSize()] - motif->getPosSetSize()*(length-1);
	if(Global::revcomp){
		maxSites = Global::posSet->total[motif->getPosSetSize()] - 2*motif->getPosSetSize()*(length-1);
	}
}


#endif
