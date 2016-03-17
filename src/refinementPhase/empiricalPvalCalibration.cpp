#include "empiricalPvalCalibration.h"
#include "StartPosUtils.h"

empiricalPvalCalibration::empiricalPvalCalibration(){
	scoreList 			= (sorted_sites*)malloc(Global::negSet->total[Global::negSet->nent] / 5 * sizeof(double));
	for(int i=0; i<Global::negSet->total[Global::negSet->nent] / 5; i++){
		scoreList[i]    = (sorted_sites)malloc(sizeof(sorted_sites_type));
	}
}

empiricalPvalCalibration::~empiricalPvalCalibration(){
	for(int i=0; i<Global::negSet->total[Global::negSet->nent]/5; i++){
		free(scoreList[i]);
	}
	free(scoreList);
}


void empiricalPvalCalibration::updateRefScoreList(float smax){
	int motifLength = M->getMotifLength();
	int firstMotifColumn = M->getFirstMotifColumn();
	double** pwm = M->getPWM();

	/* fill array scores with motif scores for each test sequence */
	refSetNb = 0;
	for(int32_t i=1; i <= Global::negSet->nent; i++){
		uint8_t* S = Global::negSet->entity[i]->S[0];
		int width = Global::negSet->entity[i]->n - motifLength + 1;
		for(int32_t j=1; j <= width; j++){
			double score = 0;
			for(int pos=0; pos < motifLength; pos++){
				uint8_t base = S[j+pos];
				if(base == 0){
					score = std::numeric_limits<double>::infinity();
					break;
				}
				score -= pwm[pos+firstMotifColumn][base];
			}

			/*if score would be -inf go to the next motif position */
			if(score > smax){ continue; }

			//for(int pos=0; pos<motifLength; pos++){
			//	fprintf(stderr, "%c", AlphaChar(S[j+pos], Global::A));
			//}
			//fprintf(stderr, "\t%s %d/%d: %f\n", Global::negSet->entity[i]->info[0], i, j, score);

			if(updateMultOcc(scoreList, refSetNb, i, j, motifLength, score, Global::negSet)){
				refSetNb++;
			}
		}
	}
	sortSitesByScore(scoreList, refSetNb);

	//fprintf(stderr, "best: %f\n", scoreList[0]->score);
	//fprintf(stderr, "worst (%d): %f\n", nbRef, scoreList[nbRef-1]->score);
	//fprintf(stderr, "\nmax score: %f\n", smax);
}

/* calculate P-value with order statistics for a specific number of motifs with a score list of
 * negative reference motifs */
double empiricalPvalCalibration::checkMotifOnRefScoreList(float smax){
	int listLength = Global::negSet->total[Global::negSet->nent] - Global::negSet->nent*M->getMotifLength();
	int k = 0;

	// best element in list is better or equal than score
	//fprintf(stderr, "\n\nsmax: %f >= scoreList[0]: %f\n", smax, scoreList[0]->score);
	if(scoreList[0]->score <= smax){
		int stepsize = 4000;
		int current = 0;
		while(k==0){
			while(current+stepsize > refSetNb){
				stepsize /= 2;
			}
			//fprintf(stderr, "current+stepsize: %d, score: %f\n", current+stepsize, scoreList[current+stepsize]->score);
			if(current+stepsize != refSetNb && scoreList[current+stepsize]->score <= smax){
				current += stepsize;
				if(current == refSetNb){
					k= current;
				}
			}else if(stepsize > 1){
				stepsize/= 2;
			}else{
				k=std::min(refSetNb, current+1);
			}
		}
	}
	//fprintf(stderr, "\tcounts: %d\t", k);
	if(k < 10) return -1;

	k = listLength-k;
	return 1 - (k+0.5)/(listLength+1); // probability of choosing a random sequence better than the one tested
}
