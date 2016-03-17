#include "LogTable.h"
#include "Globals.h"

double*   	LogTable::LOG_i = NULL;				    /* precalculated logs */
double*	  	LogTable::LOG_sum_i = NULL;			    /* precalculated sum of log(1) + .. + log(i) */
double 	  	LogTable::LOG1_1000;				    /* precalculated log(1.0/1000) */
double 	 	LogTable::LOG_Neff_pwm;						/* precalculatd log(G->neff) */
double 		LogTable::LOG_10;			    		/* precalculated log(10) */
double 		LogTable::LOG_2;			    		/* precalculated log(2) */
double 		LogTable::LOG_LOG_Seq;					/* precalculated log(log(nSeq(S))) */

double* 		LogTable::motifNbCorrection = NULL;

LogTable::LogTable(){
	LOG1_1000 = log(1.0/1000);    /* precalculated log(1.0/1000) */
	LOG_Neff_pwm = log(Global::neff_pwm);	/* precalculatd log(G->neff) */
	LOG_10 = log(10);
	LOG_2 = log(2);
	LOG_LOG_Seq = log(log(Global::posSet->nent));	/* precalculated log(log(nSeq(S))) */


	int size = std::max(10000, Global::posSet->total[Global::posSet->nent]+1);
	LOG_i = (double*)malloc(size * sizeof(double));
	LOG_sum_i = (double*)malloc(size * sizeof(double));

	LOG_i[0] = log(1e-100);
	LOG_sum_i[0] = 0;
	for(int i=1; i< size; i++){
		LOG_i[i] = log(i);
		LOG_sum_i[i] = LOG_i[i] + LOG_sum_i[i-1];
	}
}

LogTable::~LogTable(){
	free(LOG_i);
	free(LOG_sum_i);
	if(motifNbCorrection != NULL)free(motifNbCorrection);
}

void LogTable::setMotifNbCorrection(){
	motifNbCorrection = (double*)calloc(1000, sizeof(double));

	for(int K = 2; K<=1e9; K*=2){
		motifNbCorrection[K] = 0;

		for(int a=1; a<=4; a++){
			double lambda = K * Global::negBg[a];
			double beta = sqrt(K) * Global::negBg[a];
			//fprintf(stderr, "K: %d, lambda: %f, beta: %f\n", K, lambda, beta);

			double quotient = 1;
			double sum = quotient * log( (1+beta) / (lambda+beta) );
			//fprintf(stderr, "\tk: %d, quotient: %e, sum: %e\n", 0, quotient, sum);

			for(int k=1; quotient > 1e-20; k++){
				quotient *= ( lambda / k );
				sum += quotient * log( (k+1+beta) / (lambda+beta) );
				//fprintf(stderr, "\tk: %d, quotient: %e, sum: %e\n", k, quotient, sum);
			}
			motifNbCorrection[K] += lambda * exp(-lambda) * sum;
		}
		fprintf(stderr, "%d\t%f (%f)\n", K, motifNbCorrection[K], exp(motifNbCorrection[K]));
	}
	fprintf(stderr, "exit in setMotifNbCorrection\n");
	exit(-1);
}
