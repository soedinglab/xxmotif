#include "Compare_Motifs.h"

bool MotifComparison::DEBUG = false;

/***************************+
 * decide whether two motifs are similar and for what merging distance
 **************************/
int MotifComparison::calcBestMergingDist(Motif* model1, Motif* model2, int &bestRevcomp){
	/***
	 * check whether intermotif distance is below threshold
	 ***/

	bool debug = false;
	//if((model1->isTracked() || model2->isTracked()) && model1->getTotalSites() == 14 && model2->getTotalSites() == 7){
	//	cerr << *model1;
	//	cerr << *model2;
	//	debug = true;
	//}

	bool found = false;
	const int maxOffset = 10;
	double distArray[2][2*maxOffset+1];
	const double max_PWM_diversity = 0.25;

	static double probPWM1[31][21]; // max motif size
	static double probPWM2[31][21];
	// copy pwms and convert to probabilites
	for(int i=1; i<=PWM_LENGTH; i++){
		for(int j=1; j<=nAlpha(Global::A); j++)	probPWM1[i][j] = exp(model1->getPWM()[i][j]);
	}
	for(int i=model2->getFirstMotifColumn(); i <= model2->getLastMotifColumn(); i++){
		for(int j=1; j<=nAlpha(Global::A); j++) probPWM2[i][j] = exp(model2->getPWM()[i][j]);
	}

	int revComp = 0;
	if(Global::revcomp){ revComp = 1; }
	for(int revc = 0; revc <= revComp; revc++){
	  // check all possible merging distances
	  for(int offset = -maxOffset; offset <= maxOffset; offset++){
		//if(debug && revc == 0 && offset == 0) DEBUG = true;
		//else DEBUG = false;
	    double dist = calc_intermotif_dist(probPWM1, probPWM2, model1, model2, revc, offset, max_PWM_diversity);
	    //if(debug)fprintf(stderr, "\nrevcomp: %d\toffset: %d\tdist: %f\n", revc, offset, dist);
	    if(dist < max_PWM_diversity){
	    	distArray[revc][offset+maxOffset] = dist;
	    	found = true;
	    }else{
	    	distArray[revc][offset+maxOffset] = -1;
	    }
	  }
	}

	if(debug){
		fprintf(stderr, "found: %d\n", found);
		exit(-1);
	}

	if(!found) return -10000;

	double identical_minMotif_threshold = 0.2;
	if(Global::mergeMode == HIGH){
		identical_minMotif_threshold = 0.4;
	}
	//const double identical_maxMotif_threshold = 0.2;

	const int tolerance = 0;
	int bestMerge = -1;
	int best_total_identicalMotifs = 0;
	int sum_identicalMotifs = 0;
	for(int revc = 0; revc <= revComp; revc++){
		for(int merge = -maxOffset; merge <= maxOffset; merge++){
			if(distArray[revc][merge+maxOffset] == -1)continue;

			Motif mergedModel(*model1);

			int total_identicalMotifs;
			if(revc == 0){
				total_identicalMotifs = fillMergedModel(model1, model2, &mergedModel, merge, tolerance);
			}else{
				total_identicalMotifs = fillMergedModelRevcomp(model1, model2, &mergedModel, merge, tolerance);
			}
			sum_identicalMotifs += total_identicalMotifs;

			if((total_identicalMotifs > best_total_identicalMotifs) ||
			   (total_identicalMotifs == best_total_identicalMotifs && abs(merge) < abs(bestMerge)) ){
				best_total_identicalMotifs = total_identicalMotifs;
				bestMerge = merge+maxOffset;
				bestRevcomp = revc;
			}
			if(debug){
				cerr << "mergedModel: " << mergedModel;
			}
		}
	}

	if(debug){
		fprintf(stderr, "\nmerge: %d/%d bestDist: %f, identicalMotifs: %f%%\n", bestMerge, bestRevcomp, distArray[bestRevcomp][bestMerge], best_total_identicalMotifs*100.0/std::min(model1->getTotalSites(), model2->getTotalSites()));
		fprintf(stderr, "identic: %d, motifs1: %d, motifs2: %d\n", best_total_identicalMotifs, model1->getTotalSites(), model2->getTotalSites());
	}


	if(Global::mergeMode == LOW){
		if(sum_identicalMotifs < identical_minMotif_threshold * model2->getTotalSites()){
			return -10000;
		}
	}else{
		if(sum_identicalMotifs < identical_minMotif_threshold * std::max(model1->getTotalSites(), model2->getTotalSites())){
			return -10000;
		}
	}
	//if(model2->isTracked() || model1->isTracked()){
	//	fprintf(stderr, "\nmerge: bestDist: %f, identicalMotifs: %f%%\n", distArray[bestRevcomp][bestMerge], best_total_identicalMotifs*100.0/std::min(model1->getTotalSites(), model2->getTotalSites()));
	//	fprintf(stderr, "identic: %d, motifs1: %d, motifs2: %d\n", best_total_identicalMotifs, model1->getTotalSites(), model2->getTotalSites());
	//}

	return bestMerge;
}



/*************************
 * merge PWMs
 *************************/
void MotifComparison::mergePWMwithPWM(Motif* &model1, Motif* model2, int merge, int mergeRevcomp){
	/* create model for merged Motifs as copy of motif 1 */
	const int minOffset = -10;

	Motif* mergedModel = new Motif(*model1);

	if(mergeRevcomp == 0){
		fillMergedModel(model1, model2, mergedModel, merge + minOffset, 0);
	}else{
		fillMergedModelRevcomp(model1, model2, mergedModel, merge + minOffset, 0);
	}

	if(model1->isTracked() || model2->isTracked()){
		mergedModel->setTracked();
		//cerr << *model1 << *model2;
		//cerr << "=== merged Model: ===\n";
		//cerr << "merge: " << merge << "/" << mergeRevcomp << endl;
		//cerr << *mergedModel;
	}


	/* filter merged motifs that overlap with the sequence ends or that are not within posSet sequences*/
	mergedModel->filter_fitToSequence();

	/* filter merged motifs that are not within the enriched region */
	if(Global::usePositionalProbs){
		mergedModel->updateEnrichedRegion();
		//mergedModel->filter_region();
	}

	if(Global::multipleOccurrence){
		mergedModel->updatePWM_OOL_initial(Global::pseudo);
	}else{
		mergedModel->filter_oops_updatePWM(); /* filter to only one occurrence per sequence and build up PWM */
	}

	if(mergedModel->getPosSetSize() > 0){
		/* calculate pValue for merged model */
		mergedModel->updatePval(false);
		if(Global::multipleOccurrence){
			mergedModel->updatePWM_OOL_initial(Global::pseudo);
		}else{
			mergedModel->filter_oops_updatePWM(); /* filter to only one occurrence per sequence and build up PWM */
		}

		/* calculate pValue for merged model */
		mergedModel->updatePval(false);

		/* if merged model worse than the not merged ones, do not merge */
		if(MotifComparison::DEBUG){
			cerr << *mergedModel;
			fprintf(stderr, "oldPval: %f, newPval: %f\t =>", model1->getPval(),  mergedModel->getPval());
			exit(-1);
		}
	}

	//if(mergedModel->getPval() < -1000){
	//	cerr << *model1 << *model2 << *mergedModel;
	//	mergedModel->updatePval(true);
	//	exit(-1);
	//}

	if(mergedModel->getStartPosList().size() == 0 || mergedModel->getPval() >= model1->getPval()){
		if(MotifComparison::DEBUG) fprintf(stderr, "do not merge\n\n");
		if(model2->isTracked()) {
			cerr << "\n++++++Tracked Motif makes pValue worse => do not merge";
			cerr << "\n      tracked Motif switches to the other one with better pValue";
			cerr << "\n      merge: " << merge << " mergeRevcomp: " << mergeRevcomp << endl;
			cerr << "mergedModel: " << endl;
			cerr << *mergedModel;
			cerr << "model1" << endl;
			cerr << *model1;
			cerr << "      delete: " << endl;
			cerr << *model2;
			model1->setTracked();
		}
		delete mergedModel;
	}
	/* merge => substitute motif 1 with mergedModel*/
	else{
		if(MotifComparison::DEBUG) fprintf(stderr, "merge\n\n");

		if(Global::multipleOccurrence){
			mergedModel->updatePWM_OOL_initial(Global::pseudo);
		}else{
			/* filter to only one occurrence per sequence and build up PWM */
			mergedModel->filter_oops_updatePWM();
		}
		if(mergedModel->isTracked()){
			cerr << "Merge Model => " << *mergedModel;
			cerr << "delete:" << endl;
			cerr << *model1;
			cerr << "delete:" << endl;
			cerr << *model2;
		}
		delete model1;
		model1 = mergedModel;
	}
	//cerr << "done" << endl;
}

void MotifComparison::removeOverlappingMotifOccurrencies(Motif* model1, Motif* model2, int merge){
	StartPosContainer& sites = model1->getStartPosList();
	StartPosContainer& sites2 = model2->getStartPosList();

 	StartPosContainer::iterator it = sites.begin();

 	double* bg_log = Global::negSet != NULL ? Global::negBg_log : Global::posBg_log;

 	int firstMotifColomn_2 = model2->getFirstMotifColumn();
 	StartPosContainer::const_iterator end = sites.end();
 	for (StartPosContainer::iterator it2 = sites2.begin(); it2!=sites2.end(); it2++){
 		int32_t seq_2 = it2->seq;
 		int32_t pos_2 = it2->pos + firstMotifColomn_2;

	 	/* check whether motif already in list */
		while(it != end && it->seq < seq_2){it++;}
		while(it != end && it->seq == seq_2 && it->pos < pos_2-merge){it++;}
		if(it == end){	continue; }
		if(it->seq == seq_2){
			if(pos_2 - it->pos == merge){
				//cerr << "found identical motif" << endl;
				//int offset = Global::posSet->entity[it->seq]->n - Global::downstream;
				//cerr << "seq: " << it->seq << " startPos: " << it->pos + model1->getFirstMotifColumn() -1 - offset<< endl;

				uint8_t* S = Global::posSet->entity[it->seq]->S[0];
				double** pwm1 = model1->getPWM();
				double score1 = 0;
				for(motif_columns_type::const_iterator it_column = model1->getMotifColumns().begin();
						it_column != model1->getMotifColumns().end(); it_column++){
					int pos = it->pos + *it_column - 1;
					int base = S[pos];
					score1 += pwm1[*it_column][base] - bg_log[base];
					//cerr << AlphaChar(base, Global::A);
				}
				score1 /= static_cast<int>(model1->getMotifColumns().size());
				//cerr << "\tscore1: " << score1 << endl;
				double** pwm2 = model2->getPWM();
				double score2 = 0;
				for(motif_columns_type::const_iterator it_column = model2->getMotifColumns().begin();
						it_column != model2->getMotifColumns().end(); it_column++){
					int pos = it2->pos + *it_column - 1;
					int base = S[pos];
					score2 += pwm2[*it_column][base] - bg_log[base];
					//cerr << AlphaChar(base, Global::A);
				}
				score2 /= static_cast<int>(model2->getMotifColumns().size());
				//cerr << "\tscore2: " << score2 << endl;

				/* motif fits better to second pwm */
				if(score2 > score1){
					it = sites.erase(it);
				/* motif fits better to first pwm */
				}else{
					it2 = sites2.erase(it2);
				}
			}
		}
	}
 	model1->updatePWM_OOL_initial(Global::pseudo);
 	model2->updatePWM_OOL_initial(Global::pseudo);
 }




