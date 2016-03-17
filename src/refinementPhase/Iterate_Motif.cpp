#include "Iterate_Motif.h"

IterateMotif::IterateMotif(){
	elongPWM = (double**)malloc(31*sizeof(double*));
	bestElongPWM = (double**)malloc(31*sizeof(double*));
	for(int i=0; i<= PWM_LENGTH; i++){
		elongPWM[i] = (double*)calloc(nAlpha(Global::A)+1, sizeof(double));
		bestElongPWM[i] = (double*)calloc(nAlpha(Global::A)+1, sizeof(double));
	}

	skipRoundThreshold = log(0.7);
}

IterateMotif::~IterateMotif(){
	for(int i=0; i<= PWM_LENGTH; i++){
		free(elongPWM[i]);
		free(bestElongPWM[i]);
	}
	free(elongPWM);
	free(bestElongPWM);
}

void IterateMotif::optimizeMotif(Motif* motif, bool keepPosSetSize, bool plus){
	int bestMotifNb; double bestPval;

	for(int iter = 0; iter < 2; iter++){
			Kmer_Generator_Reloaded<dna_hash_t> &kmerHash = Kmer_Generator_Reloaded<dna_hash_t>::getInstance();
			kmerHash.reinitialize(motif->getPWM(), motif->getMotifColumns(), motif->getPosSetSize(), Global::posSet->total[motif->getPosSetSize()]);
			StartPosUpdater::getInstance().fill_sites_with_possible_startPositions(kmerHash, keepPosSetSize, motif);
		//cout << *motif;
		StartPosUpdater::getInstance().update_PWM_with_bestSubset(motif, bestMotifNb, bestPval, plus, Global::pseudo);
	}

	StartPosUpdater::getInstance().setPositionalPval(motif);
	//StartPosUpdater::getInstance().update_pwm_with_bestSubset_All(motif, pseudo);
	//fprintf(stderr, "bestMotifNb, %d, bestPval: %e\n", bestMotifNb, bestPval);
}


bool IterateMotif::maximizeMotifLength(Motif* motif, bool lengthCorrection, bool palin){
	//motif->printFullPWM(cerr);

	/* motif does not have more than maxMatchPositions match positions */
	int bestMotifNb = 0;
	bool noChange = false;
	double bestPval = std::numeric_limits<double>::max();
	int firstMotifColumn = motif->getFirstMotifColumn();
	StartPosContainer bestStartPosList, realStartPosList;

	/* try each change of match column */
	int frontStart = lengthCorrection ? -2 : 0;
	//int frontEnd=0;int backStart=0; int backEnd=0;
	int frontEnd   = lengthCorrection ?  2 : 0;
	int backStart  = lengthCorrection ? -2 : 0;
	int backEnd    = lengthCorrection ?  2 : 0;

	motif_columns_type bestElongation = motif->getMotifColumns();
	for(int front = frontStart; front <= frontEnd; front ++){
		for(int back = backStart; back <= backEnd; back++){
			//fprintf(stderr, "front: %d, back: %d\n", front, back);
			motif_columns_type tmpColumns;
			bool skipRound;
				skipRound = fillTmpColumn_noGaps(front, back, tmpColumns, motif);
			if(skipRound){continue; }

			//cerr << "front: " << front << " back: " << back << " : ";
			//for(motif_columns_type::iterator it_column = tmpColumns.begin();it_column != tmpColumns.end(); it_column++)
			//	cerr << *it_column << " ";
			//cerr << endl << "before" << endl;
			//printPWM(motif->getPWM(), tmpColumns, cerr);

			/* calculate pVal for new match column combination */
			int motifNb=0;
			double overrepPval=1e100;

			for(int iter = 0; iter < 3; iter++){
				double **pwm = elongPWM;
				if(iter == 0)pwm = motif->getPWM();

				/* calculate pValues for optimized matrix */
					Kmer_Generator_Reloaded<dna_hash_t> &kmerHash = Kmer_Generator_Reloaded<dna_hash_t>::getInstance();
					kmerHash.reinitialize(pwm, tmpColumns, motif->getPosSetSize(), Global::posSet->total[motif->getPosSetSize()]);
					StartPosUpdater::getInstance().fill_elongSites_with_possible_startPositions(kmerHash, motif, pwm, tmpColumns, palin);

				bool plus     = (iter == 0) ?  true : false;
				StartPosUpdater::getInstance().update_elongPWM_with_bestSubset(motif, elongPWM, tmpColumns, motifNb, overrepPval, plus, Global::pseudo);

				if(motifNb == 0)break;
				//cerr << endl;
				//if(iter!=0)fprintf(stderr, "overrepresentation Pval: exp(%e) = %e (%d)\n", overrepPval, exp(overrepPval), motifNb);
				//else fprintf(stderr, "motifNbPlus: %d\n", motifNb);
				//printPWM(elongPWM, tmpColumns, cerr);
			}
			if(motifNb == 0)continue;

			//double newPval = StartPosUpdater::getInstance().getCombinedPval(motif, elongPWM, tmpColumns, motifNb, overrepPval);
			double LOG_Bonferonni = calculate_log_bonferonni(tmpColumns, LogTable::LOG_Neff_pwm);
			double newPval = overrepPval + LOG_Bonferonni;

			//cerr << endl;
			//printPWM(elongPWM, tmpColumns, cerr);
			//fprintf(stderr, "posSetSize: %d, best combined: %e (%d)\n\n", motif->getPosSetSize(), exp(newPval), motifNb);
			//exit(-1);
			if(newPval < bestPval){
				bestPval = newPval;
				bestElongation = tmpColumns;
				bestMotifNb = motifNb;
				if(front == 0 && back == 0) noChange = true;
				else noChange = false;

				/* copy current elong pwm to best elong pwm */
				for(motif_columns_type::iterator it_col = tmpColumns.begin(); it_col != tmpColumns.end(); it_col++){
					for(int j = 1; j<= nAlpha(Global::A); j++){
						bestElongPWM[*it_col][j] = elongPWM[*it_col][j];
					}
				}
				motif->setWorstScore(StartPosUpdater::getInstance().fill_startPosList_with_elongSites(realStartPosList, bestMotifNb, tmpColumns));
			}
			//fprintf(stderr, "done\n");
		}
	}
	if(bestMotifNb == 0) return noChange;

	motif->setPval(bestPval);
	motif->setMotifColumns(bestElongation);

	motif->getStartPosList().clear();
	firstMotifColumn = (motif->getLength() - motif->getMotifLength()) / 2 + 1;
	int offset = firstMotifColumn - motif->getFirstMotifColumn();

	int count = 0;
	for(StartPosContainer::iterator it = realStartPosList.begin(); it != realStartPosList.end(); it++, count++){
		it->pos = static_cast<int32_t>(it->pos - offset);
		motif->getStartPosList().push_back(*it);
	}

	/* copy best elong pwm to motif pwm */
	for(motif_columns_type::iterator it_col = bestElongation.begin(); it_col != bestElongation.end(); it_col++){
		for(int j = 1; j<= nAlpha(Global::A); j++){
			motif->getPWM()[*it_col][j] = bestElongPWM[*it_col][j];
		}
	}

	if(offset != 0)	motif->offset(offset);

	StartPosUpdater::getInstance().setPositionalPval(motif);

	return noChange;
}


void IterateMotif::printPWM(double** pwm, motif_columns_type& motifColumn,  ostream &os) const{
	os << endl << "PWM" << endl << "\t\t";
	for(int i=0;i<motifColumn.back()-motifColumn.front()+1;i++)
		os << std::setprecision(2) << std::setw(6) << i+1 << "\t";
	os << endl << "\t-------";
	for(int i=motifColumn.front();i<=motifColumn.back();i++)os << "--------";
	for(int i=1;i<=nAlpha(Global::A);i++){
		os << endl << "\t" << AlphaChar(i,Global::A) << "\t";
		os << std::resetiosflags( std::ios::scientific );
		os << std::fixed;
		motif_columns_type::const_iterator it = motifColumn.begin();
		while(true){
			int column = *it;
			os << std::setw(6) << (int)(exp(pwm[*it][i]) * 10000)/100.00 << "\t";
			it++;
			if(it == motifColumn.end())break;
			for(int col = column+1; col < *it; col++) os << std::setw(6) << "-" << "\t";
		}
	}
	os << endl << endl;
}

