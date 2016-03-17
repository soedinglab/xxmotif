#include "backgroundDistribution.h"

#include <algorithm>
#include <cassert>
#include "utils.h"
#include "NullModel.h"

#include "em/hoNullModel.h"

void setBackgroundDistribution(){
	if( Global::negSet != NULL ){
		NullModel::init( Global::negSet, Global::pseudocountsFactor, Global::countsOffset, Global::order, true, NULL );
	} else{
		NullModel::init( Global::posSet, Global::pseudocountsFactor, Global::countsOffset, Global::order, true, NULL );
	}

	if( Global::em ){

		ss_type sequences;
		sequences = Global::negSet == NULL ? Global::posSet : Global::negSet;

		hoNullModel::init( sequences, Global::alphaBg, Global::modelOrderBg, Global::freqs );
		if( Global::saveModels ){
			fprintf( stderr, "%s\n", baseFileName( sequences->name ) );
			hoNullModel::save( baseFileName( sequences->name ) );
		}
	}
}

void checkSequenceSet(){
	int len, base, base2, base3;
	double totPos=0;
	double totNeg = 0;
	const int aS = nAlpha(Global::A) + 1;
	int ***trimersPosSet = new int**[aS];
	int ***trimersNegSet = new int**[aS];
	for(int i=0; i<aS; i++){
	trimersPosSet[i] = new int*[aS];	
	trimersNegSet[i] = new int*[aS];
	for(int j=0; j<aS; j++){
	    trimersPosSet[i][j] = new int[aS];
	    trimersNegSet[i][j] = new int[aS];
		for(int k=0; k<aS; k++){
			trimersPosSet[i][j][k] = 0;
			trimersNegSet[i][j][k] = 0;
		}
	    }
	}

	if(Global::negSet != NULL){
		for(int i=1; i <= Global::negSet->nent; i++){
			len = Global::negSet->entity[i]->n;
			for(int k=0; k<1; k++){
				totNeg += len;
				base = (int)Global::negSet->entity[i]->S[k][1];
				base2= (int)Global::negSet->entity[i]->S[k][2];
				Global::negBg[base]++;
				Global::negBg[base2]++;
				for(int j=3; j <= len; j++){
					base3= (int)Global::negSet->entity[i]->S[k][j];					
					Global::negBg[base3]++;
					trimersNegSet[base][base2][base3]++;
					base=base2;
					base2=base3;
				}
			}
		}
	}	

	for(int i=1; i <= Global::posSet->nent; i++){
		len = Global::posSet->entity[i]->n;
		for(int k=0; k<1; k++){
			totPos += len;
			base = (int)Global::posSet->entity[i]->S[k][1];
			base2= (int)Global::posSet->entity[i]->S[k][2];
			Global::posBg[base]++;
			Global::posBg[base2]++;
			
			for(int j=3; j <= len; j++){
				base3 = (int)Global::posSet->entity[i]->S[k][j];
				Global::posBg[base3]++;
				trimersPosSet[base][base2][base3]++;
				base=base2;
				base2=base3;
			}
		}
	}
	for(int i=1;i<aS;i++){
		if(Global::negSet != NULL){
			Global::negBg[i] = std::max(1e-100, Global::negBg[i] / (totNeg-Global::negBg[0]));
			Global::negBg_log[i] = log(Global::negBg[i]);
		}
		Global::posBg[i] = std::max(1e-100, Global::posBg[i] / (totPos-Global::posBg[0]));
		Global::posBg_log[i] = log(Global::posBg[i]);
	}

	printf("\nCheck composition in sequence sets:\n");
	printf("posSet:\t");
	double sum = 0;
	for(int i=1;i<aS;i++){
		printf("%c: %.2f%%\t", AlphaChar(i, Global::A), Global::posBg[i]*100);
		sum += Global::posBg[i];
	}
	printf("\n");
	//printf("Σ: %.2f%%\n", sum * 100);
	if(Global::negSet != NULL){
		printf("%s:\t", "negSet");
		double sum = 0;
		for(int i=1;i<aS;i++){
			printf("%c: %.2f%%\t", AlphaChar(i, Global::A), exp(Global::negBg_log[i])*100);
			sum += exp(Global::negBg_log[i]);
		}
		//printf("Σ: %.2f%%\n", sum * 100);
		printf("\n");
		printf("\n");

		/* caluclate RMSD of timer counts without N's */
		int posCount = 0;
		int negCount = 0;
		for(int i=1; i<aS; i++){
			for(int j=1; j<aS; j++){
				for(int k=1; k<aS; k++){
					//printf("%c%c%c: %d, %d\n", AlphaChar(i, Global::A), AlphaChar(j, Global::A), AlphaChar(k, Global::A), trimersPosSet[i][j][k], trimersNegSet[i][j][k]);
					posCount += trimersPosSet[i][j][k];
					negCount += trimersNegSet[i][j][k];
				}
			}
		}
		//printf("\n\nposCount : %d , negCount: %d \n", posCount, negCount);

		double sumSquaredDists = 0;
		int N = 0;
		for(int i=1; i<aS; i++){
			for(int j=1; j<aS; j++){
				for(int k=1; k<aS; k++){
					//double dist = (trimersPosSet[i][j][k]*100.0)/posCount - (trimersNegSet[i][j][k]*100.0)/negCount;
					double dist = (trimersPosSet[i][j][k]*1.0)/posCount - (trimersNegSet[i][j][k]*1.0)/negCount;
					//printf("%c%c%c: %f, %f\n", AlphaChar(i, Global::A), AlphaChar(j, Global::A), AlphaChar(k, Global::A), (trimersPosSet[i][j][k]*1.0)/posCount, (trimersNegSet[i][j][k]*1.0)/negCount);					
					sumSquaredDists += dist*dist;
					N++;
				}
			}
		}
 		//printf("sumSquared: %f, N: %d\n", sumSquaredDists, N);
		double RMSD = sumSquaredDists / N;
		RMSD = sqrt(RMSD);
		//printf("Root mean square deviation of trimers in positive and negative set: %f\n\n", RMSD);
		/* Evalutaed on Harbison data */
		double mean = 0.001917581;
		double sd = 0.0009103165;	
		if(RMSD > mean + 2*sd){
			printf("  =====================================================================\n");
			printf("  ============================== WARNING ==============================\n");
			printf("  ==  The RMSD between the trimer probabilities within the positive  ==\n");
			printf("  ==      and negative set is %4.1f sd higher than recommended!       ==\n", (RMSD-mean)/sd);
			printf("  ==   --> The results might be better without using a negative set  ==\n");
			printf("  ==  (Bad negative sets lead to long run times and motifs that are  ==\n");
			printf("  ==                  too long and too significant!)                 ==\n");
			printf("  =====================================================================\n");
			printf("  =====================================================================\n\n");
		}
	}

  for(int i=0; i<aS; i++){
    for(int j=0; j<aS; j++){
      delete[] trimersPosSet[i][j];
      delete[] trimersNegSet[i][j];
    }
    delete[] trimersPosSet[i];
    delete[] trimersNegSet[i];
  }
  delete[] trimersPosSet;
  delete[] trimersNegSet;

	//exit(-1);
}

