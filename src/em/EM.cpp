#include "EM.h"
#include "hoNullModel.h"
#include "../NullModel.h"
#include "../output.h"
#include "../refinementPhase/Motif.h"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;

void EM::go(){

	bool printModels = false;
	bool printPlacements = false;
	bool printPosteriors = false;

	bool checkLikelihoods = false;

	int N = Global::posSet->nent;

	if( Global::verbose ){
		printf( " ____\n"
				"|    |\n"
				"| EM |\n"
				"|____|\n"
				"\n"
				"on %d sequences\n", N );
	}

	int b, i, k, l, n, w, x;
	int LW1, LW2, W;

	int Worder;

	// count iterations
	int round;
	// convergence flag
	bool iterate;
	// likelihood or model parameter difference between EM iterations
	double difference;

	// posterior probabilities
	double **counts;
	double **countsx;
	// conditional probabilities
	double **conds;
	// probabilities
	double **probs;
	// l(ast) iteration's probabilities
	double **lprobs;

	// background probabilities over sequence positions
	double** bg;
	// posterior probabilities over sequence positions
	double** posterior;
	// likelihood ratios over sequence positions
	double** likelihoods;
	// model placements without N character
	bool** placements;
	// background model shifts
	// due to N character before model placements
	int** shifts;
	int shift;

	// parameters for optimizing alpha
	double Qfunc_all;
	double Qfunc_grad;
	double Grad_all;

	// likelihood ratio at state k
	double likelihood_ratio;
	// likelihood at state k
	double likelihood_k;
	// likelihood
	double likelihood;
	// l(ast) iteration's likelihood
	double llikelihood;
	// normalization factor
	double normFactor;

	int order;
	float q;
	int fields_maxorder;

	std::stringstream str;
	FILE* f_likelihoods = NULL;


	// to calculate sequence-specific likelihoods
	double likelihood_sequence;
	// to calculate estimated complete log-likelihood Q + log prior
	double bound;
	// l(ast) iteration's bound
	double lbound;
	// to calculate sequence-specific bound
	double bound_sequence;

	double prior_state;
	double prior_state_0;

	ss_type sequences = Global::posSet;
	unsigned char* sequence;
	float weight;

	int L = -1; // sequence lengths must not differ
	for( int n=1; n <= N; n++ ){

		// ignore remaining sequences in S if any
		sequence = sequences->entity[n]->S[0];

		if( L < 0 ){
			L = sequences->entity[n]->n;
		} else if( L != sequences->entity[n]->n ){
			fprintf( stderr, "Sequence lengths must not differ.\n" );
			exit(-1);
		}
	}

	int modelCounter = 1;

	Motif* motif;
	list<Motif*>::const_iterator iter;
	for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end();
			iter++ ){

		if( Global::verbose ){
			printf( " _______\n"
					"|       |\n"
					"| MODEL |\n"
					"|_______|\n"
					"\n"
					"no. %d\n", modelCounter );
		}

		motif = *iter;

		order = motif->getOrder();
		q = motif->getSpecificityFactor();

		b = motif->getFirstMotifColumn();
		W = motif->getMotifLength();

		if( W > L ){
			fprintf( stderr, "Sequences shorter than model no. %d\n",
					modelCounter );
			continue;
		}

		// more readable names
		counts = motif->getCounts();
		countsx = motif->getCountsx();

		conds = motif->getConds();
		probs = motif->getPWM();

		// starting index to access max. order probabilities
		fields_maxorder = motif->getFields() - motif->getCoefficients()[order+1]
		                                                                + 1;

		// precalculations
		LW1 = L-W+1;
		LW2 = L-W+2;
		Worder = W-order;

		// (re)settings
		round = 0;
		iterate = true;
		difference = DBL_MAX;
		llikelihood = 0;
		lbound = 0;

		if( Global::revcomp ){
			prior_state = static_cast<double>( q ) / static_cast<double>(
					L-W+1-W+1 );
		} else{
			prior_state = static_cast<double>( q ) / static_cast<double>( LW1 );
		}
		prior_state_0 = 1.0 - static_cast<double>( q );

		// allocate data structures
		bg = ( double** )malloc( (N+1)*sizeof( double* ) );
		likelihoods = ( double** )malloc( (N+1)*sizeof( double* ) );
		placements = ( bool** )malloc( (N+1)*sizeof( bool* ) );
		posterior = ( double** )malloc( (N+1)*sizeof( double* ) );
		shifts = ( int** )malloc( (N+1)*sizeof( int* ) );

		for( n=1; n <= N; n++ ){

			// allocate data structures for positions
			bg[n] = ( double* )calloc( LW2, sizeof( double ) );
			likelihoods[n] = ( double* )calloc( LW2, sizeof( double ) );
			placements[n] = ( bool* )calloc( LW2, sizeof( bool ) );
			posterior[n] = ( double* )calloc( LW2, sizeof( double ) );
			shifts[n] = ( int* )calloc( LW2, sizeof( int ) );

			// precalculate model placements and background model shifts
			// across sequences
			// across sequence positions
			// for subsequences
			sequence = sequences->entity[n]->S[0];
			for( k=1; k <= LW1; k++ ){
				for( l=0; l < W; l++ ){
					if( !( sequence[k+l] ) ){
						break; // N character
					}
				}
				if( l == W ){
					// placement possible
					placements[n][k] = true;
					// calculate background model shifts
					if( Global::modelOrderBg > 0 ){
						shift = Global::modelOrderBg;
						l = k-1;
						while( ( l > 0 ) && ( ( k-l ) <= Global::modelOrderBg )
						){
							if( !( sequence[l] ) ){
								break; // N character
							}
							shift--;
							l--; // back in sequence
						}
						shifts[n][k] = shift;
					}
				}
			}
		}

		if( printPlacements ){
			printf( " ____________\n"
					"|            |\n"
					"| PLACEMENTS |\n"
					"|____________|\n\n" );
			for( n=1; n <= N; n++ ){
				for( k=1; k <= LW1; k++ ){
					printf( "%d ", placements[n][k] );
				}
				printf( "\n" );
			}
			printf( "\n" );
			printf( " ________\n"
					"|        |\n"
					"| SHIFTS |\n"
					"|________|\n\n" );
			for( n=1; n <= N; n++ ){
				for( k=1; k <= LW1; k++ ){
					printf( "%d ", shifts[n][k] );
				}
				printf( "\n" );
			}
			printf( "\n" );
		}

		// allocate data structure
		// last iteration's probabilities
		lprobs = ( double** )malloc( W*sizeof( double* ) );
		for( l=0; l < W; l++ ){ // 0:(W-1) corresponds to b:(b+W-1)
			// allocate data structure for positions
			lprobs[l] = ( double* )calloc( motif->getFields(), sizeof( double )
			);
			// copy probabilities
			std::memcpy( lprobs[l], probs[b+l], motif->getFields()*sizeof(
					double ) );
		}

		// precalculate null model probabilities
		// across sequences
		// across sequence positions
		// for subsequences across model positions

		// calculate probabilities across sequences
		for( n=1; n <= N; n++ ){
			sequence = sequences->entity[n]->S[0];

			// calculate probabilities across positions
			for( k=1; k <= LW1; k++ ){

				if( placements[n][k] ){

					likelihood_k = 1.0;
					// across model positions
					//					for( i=k; (i-k) < W; i++ ){
					//						// using conditional probabilities
					//						likelihood_k *= hoNullModel::getConds()[
					//						                hoNullModel::sub2ind( sequence,
					//								        max( 1, i-Global::modelOrderBg ),
					//								        min( i-1, Global::modelOrderBg ) )];
					//					}
					shift = shifts[n][k];
					for( i=k; (i-k) < W; i++ ){
						// using conditional probabilities
						likelihood_k *= hoNullModel::getConds()[
						                                        hoNullModel::sub2ind( sequence,
						                                        		i-Global::modelOrderBg+shift,
						                                        		Global::modelOrderBg-shift )];
						shift = max( 0, shift-1 );
					}
					bg[n][k] = likelihood_k;
				}
			}
		}

		do{
			/**********
			 * E step: calculate posterior
			 **********/

			// calculate posterior across sequences
			for( n=1; n <= N; n++ ){
				sequence = sequences->entity[n]->S[0];

				// reset normalization factor
				normFactor = 0;

				// calculate posterior across positions
				for( k=1; k <= LW1; k++ ){

					if( placements[n][k] ){

						// using max. order probabilities
						likelihood_k = probs[b+order][motif->sub2ind( sequence,
								k, order )];
						// across remaining sequence
						for( i=k+1; (i-k) < Worder; i++ ){
							// using conditional probabilities
							likelihood_k *= conds[b+(i-k)+order][motif->sub2ind(
									sequence, i, order )];
						}

						// likelihood ratio
						likelihood_ratio = likelihood_k / bg[n][k];
						likelihoods[n][k] = likelihood_ratio;

						// likelihood ratio times prior
						posterior[n][k] = likelihood_ratio * prior_state;
						normFactor += posterior[n][k];
					}
				}

				// k=0
				likelihoods[n][0] = 1.0;
				posterior[n][0] = prior_state_0;
				normFactor += posterior[n][0];

				// posterior probability normalizations
				for( k=0; k <= LW1; k++ ){
					posterior[n][k] /=normFactor;
				}
			}

			if( printPosteriors ){
				printf( " ____________\n"
						"|            |\n"
						"| POSTERIORS |\n"
						"|____________|\n\n" );
				for( n=1; n <= N; n++ ){
					for( k=0; k <= LW1; k++ ){
						printf( "%.4f ", posterior[n][k] );
					}
					printf( "\n" );
				}
				printf( "\n" );
			}

			if( Global::saveExpectationMaximizationLikelihoods ){

				str.str( "" );
				str << Global::outputDirectory << "/" << Global::shortFileName
						<< "-" << modelCounter << "-sequence-log-likelihood-iter-"
						<< round+1 << ".txt";
				f_likelihoods = fopen( str.str().c_str(), "w" );
			}

			// total log likelihood
			likelihood = 0.0;
			for( n=1; n <= N; n++ ){

				// sequence likelihood ratio
				likelihood_sequence = 0.0;
				for( k=1; k <= LW1; k++ ){
					if( placements[n][k] ){
						likelihood_sequence += likelihoods[n][k];
					}
				}

				// likelihood ratio times prior
				likelihood_sequence *= prior_state;

				// k=0
				likelihood_sequence += likelihoods[n][0] * prior_state_0;

				// sequence log likelihood
				weight = sequences->entity[n]->weight;
				likelihood_sequence = weight * log( likelihood_sequence );
				likelihood += likelihood_sequence;

				if( Global::saveExpectationMaximizationLikelihoods ){
					fprintf( f_likelihoods, "%g\n", likelihood_sequence );
				}
			}

			if( Global::saveExpectationMaximizationLikelihoods ){

				fclose( f_likelihoods );

				str.str( "" );
				str << Global::outputDirectory << "/" << Global::shortFileName
						<< "-" << modelCounter <<
						"-sequence-state-likelihood-ratio-iter-" << round+1 <<
						".txt";
				f_likelihoods = fopen( str.str().c_str(), "w" );

				for( n=1; n <= N; n++ ){
					if( placements[n][1] ){
						fprintf( f_likelihoods, "%g", likelihoods[n][1] );
					} else{
						fprintf( f_likelihoods, "NA" );
					}
					for( k=2; k <= LW1; k++ ){
						if( placements[n][k] ){
							fprintf( f_likelihoods, " %g", likelihoods[n][k] );
						} else{
							fprintf( f_likelihoods, " NA" );
						}
					}
					fprintf( f_likelihoods, "\n" );
				}
				fclose( f_likelihoods );
			}

			if( Global::saveExpectationMaximizationModels ){
				str.str( "" );
				str << Global::shortFileName << "-" << modelCounter <<
						"-model-iter-" << round+1;
				saveInterpolatedMarkovModel( *motif, str.str().c_str() );
			}

			if( round > 0 ){

				if( Global::likelihoodConvergence ){
					// calculate log likelihood difference
					// to last iteration's likelihood
					difference = likelihood - llikelihood;
				} else{
					// calculate summed up parameter differences
					// using max. order probabilities
					difference = 0;
					for( w=0; w < W; w++ ){
						for( x=fields_maxorder; x <= motif->getFields(); x++ ){
							difference += fabs( probs[b+w][x] - lprobs[w][x] );
						}
					}
				}
			}

			if( Global::verbose ){

				printf( "Log likelihood: %10.4f", likelihood );

				if( round > 0 ){
					if( likelihood < llikelihood ){
						printf( " decreasing..." );
					}
				}
				printf( "\n" );
			}

			if( checkLikelihoods ){

				// expected complete log likelihood
				bound = 0.0;
				for( n=1; n <= N; n++ ){
					bound_sequence = 0.0;
					for( k=1; k <= LW1; k++ ){
						if( placements[n][k] ){
							bound_sequence += posterior[n][k] * ( log(
									prior_state ) + log(
											likelihoods[n][k] ) );
						}
					}
					bound_sequence += posterior[n][0] * ( log( prior_state_0 ) +
							log( likelihoods[n][0] ) );

					weight = sequences->entity[n]->weight;
					bound += weight * bound_sequence;
				}

				printf( "Expected complete log likelihood: %10.4f", bound );

				if( round > 0 ){

					if( bound < lbound ){
						fprintf( stdout, " decreasing...\n" );
					} else{
						printf( "\n" );
					}

					if( ( likelihood - llikelihood ) < ( bound - lbound ) ){
						printf( "Log likelihood difference smaller\n" );
						// compared to expected complete log likelihood
						// difference
					}
				} else{
					printf( "\n" );
				}

			}

			// check convergence
			if( difference < Global::epsilon ){

				iterate = false;
			} else{

				// last iterations's likelihood
				llikelihood = likelihood;

				// last iteration's probabilities
				for( l=0; l < W; l++ ){ // 0:(W-1) corresponds to b:(b+W-1)
					std::memcpy( lprobs[l], probs[b+l], motif->getFields()*sizeof( double ) );
				}

				if( checkLikelihoods ){
					lbound = bound;
				}

				/********
				 * M step: update model parameters
				 ********/

				// reset counts
				motif->resetCounts();
				motif->resetCountsx();

				// calculate counts across sequences
				for( n=1; n <= N; n++ ){
					sequence = sequences->entity[n]->S[0];
					weight = sequences->entity[n]->weight;

					// calculate counts across sequence positions
					for( k=1; k <= LW1; k++ ){

						if( placements[n][k] ){
							// calculate counts across model positions
							for( l=0; l < W; l++ ){
								// calculate counts across orders
								for( i=0; i <= order && l+i < W; i++ ){
									// add posterior to counts
									counts[b+l+i][
									              motif->sub2ind( sequence, k+l, i )] +=
									            		  weight * static_cast<float>( posterior[n][k]
									            		  );
									if( i > 0 ){
										// k-mer never followed by N character
										countsx[b+l+i-1][
										                 motif->sub2ind( sequence, k+l, i-1 )] +=
										                		 weight * static_cast<float>(
										                				 posterior[n][k] );
									}
								}
							}
						}
					}
				}


				if( Global::lastCondsPseudoCounts ){
					motif->setLastConds( conds );
				}

				if( Global::interpolate ){
					// update interpolated Markov model probabilities
					// keep alpha values fixed
					if( ( round+1 ) < 1 ){
						motif->calculateInterpolatedMarkovModelProbs(
								Global::lastCondsPseudoCounts );
					} else{
						motif->calculateInterpolatedMarkovModelProbs(
								Global::lastCondsPseudoCounts, 0 );
					}
				} else{
					// update Markov model probabilities
					motif->calculateMarkovModelProbs();
				}

				round++;
			}

			if( printModels ){
				printInterpolatedMarkovModel( *motif );
			}

			//printf(" round %d  ",round);
			if(round < 20){ continue;}
			//if(round == 20){
			//	printInterpolatedMarkovModel( *motif );
			//	exit(0);
			//}

			if( Global::learnHyperParameter ){

				// start only after em-algo has run for x-times


				int index [5] = {4,20,84,340,1364};
				double alpha = 0.0;
				double step_size_alpha = 1;
				int alpha_max = 1000;
				int order = 2;
				int y;
				list<int>::iterator j;

				double** condsW  = motif->getConds();
				double** countsW = motif->getCounts();
				double* ncounts = ( double* )calloc( motif->getLength()+1, sizeof( double ) );

				for( j=motif->getMotifColumns().begin(); j != motif->getMotifColumns().end(); j++ ){
					for( int i=1; i <= nAlpha( Global::A ); i++ ){
						ncounts[*j] += motif->getCounts()[*j][i];
					}
				}

				FILE * bgFile = fopen ("/home/kiesel/Projects/AwesomeProject/Eclipse_Kepler_Workspace/gimmemotif/example/probsBG.txt","w");
				for(j = motif->getMotifColumns().begin(); j != motif->getMotifColumns().end(); j++){
					for(y=1;y<=4;y++){
						fprintf(bgFile,"%.16f\t ", Global::posBg[y]);
					}
					fprintf(bgFile,"\n");
				}
				fclose(bgFile);

				char cFname[1024];
				sprintf(cFname,"/home/kiesel/Projects/AwesomeProject/Eclipse_Kepler_Workspace/gimmemotif/example/condsCheck_before_order%d.txt",order);

				FILE * cFile = fopen (cFname,"w");
				for(j = motif->getMotifColumns().begin(); j != motif->getMotifColumns().end(); j++){
					for(y=1;y<=index[order];y++){
						fprintf(cFile,"%.16f\t ",condsW[*j][y]);
					}
					fprintf(cFile,"\n");
				}
				fclose(cFile);

				char logcFname[1024];
				sprintf(logcFname,"/home/kiesel/Projects/AwesomeProject/Eclipse_Kepler_Workspace/gimmemotif/example/LogcondsCheck_before_order%d.txt",order);

				FILE * logcFile = fopen (logcFname,"w");
				for(j = motif->getMotifColumns().begin(); j != motif->getMotifColumns().end(); j++){
					for(y=1;y<=index[order];y++){
						fprintf(logcFile,"%.16f\t ",log(condsW[*j][y]));
					}
					fprintf(logcFile,"\n");
				}
				fclose(logcFile);

				char tFname[1024];
				sprintf(tFname, "/home/kiesel/Projects/AwesomeProject/Eclipse_Kepler_Workspace/gimmemotif/example/countsCheck_before_order%d.txt",order);

				FILE * tFile = fopen (tFname,"w");
				for(j = motif->getMotifColumns().begin(); j != motif->getMotifColumns().end(); j++){
					for(y=1;y<=index[order];y++){
						fprintf(tFile,"%.16f \t" ,countsW[*j][y]);
					}
					fprintf(tFile,"\n");
				}
				fclose(tFile);

				// check likelihood -> am I using the correct values? are they overwritten beforehand?
				// compare to R version
				// check range/scale from likelihood -> needs to be similar to prior!

				char pFname[1024];
				sprintf(pFname, "/home/kiesel/Projects/AwesomeProject/Eclipse_Kepler_Workspace/gimmemotif/example/Qfunc_output_order%d.txt",order);

				FILE * pFile = fopen (pFname,"w");
				fprintf(pFile,"Alpha\tQfunc_all\tQfunc_grad\tQfunc_p1\tQfunc_p2\tQfunc_p3\tQfunc_p4\tGrad_all\tGrad_p1\tGrad_p2\tGrad_p3\tzeroOrder\tpseudoOrder\n");

				for(int i = int(1/step_size_alpha); i<=alpha_max; i++ ) {
					alpha = i * step_size_alpha;

					// update conds for alpha
					uint8_t* kmer_update = ( uint8_t* )malloc( ( order+1 )*sizeof( uint8_t ) );
					motif->updateConds(kmer_update,order,alpha,0, ncounts);
					free(kmer_update);

////					if(i == 1){
//						double** condsG  = motif->getConds();
//						char gFname[1024];
//						sprintf(gFname,"/home/kiesel/Projects/AwesomeProject/Eclipse_Kepler_Workspace/gimmemotif/example/condsCheck_after_order%d_a%d.txt",order,i);
//
//						FILE * gFile = fopen (gFname,"w");
//						for(j = motif->getMotifColumns().begin(); j != motif->getMotifColumns().end(); j++){
//							for(y=1;y<=index[order];y++){
//								fprintf(gFile,"%f\t ",condsG[*j][y]);
//							}
//							fprintf(gFile,"\n");
//						}
//						fclose(gFile);
////					}

					// calculate Qfunction and Gradient for alpha
					uint8_t* kmer_score = ( uint8_t* )malloc( ( order+1 )*sizeof( uint8_t ) );
					motif->resetMotifParams(order, alpha);
					calcQandGrad(*motif, kmer_score,order,alpha,0);
					free(kmer_score);

					Qfunc_all = motif->getQfunc_p1() + motif->getQfunc_p2() + motif->getQfunc_p3() + motif->getQfunc_p4();
					Qfunc_grad = W * pow(4,order) * double(digamma(alpha + 4)) - motif->getQfunc_grad();
					Grad_all = motif->getGrad_p1() - motif->getGrad_p2() + motif->getGrad_p3();


					fprintf(pFile,"%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",alpha,Qfunc_all,Qfunc_grad,motif->getQfunc_p1(),motif->getQfunc_p2(), motif->getQfunc_p3(),motif->getQfunc_p4(), Grad_all,motif->getGrad_p1(),motif->getGrad_p2(),motif->getGrad_p3(),motif->getzeroOrder(),motif->getpseudoOrder());
				}
				fclose(pFile);
				exit(0);

				//					//explanaition for sub2ind()
				//					//first parameter  := sequence defines the beginning position within the kmer. with +1 or +2 I can reduce the kmer by "removing the first (or the forst 2) nucleotides hence reduce the kmer length
				//					//second parameter := 0, always zero since here I use the public available function , with 0 this is the same as using the private function
				//					//third parameter  := order, how many nucleotides to the right do I want to look at.
				//					// examples:
				//					// let's say kmer = "abc"
				//					// and my order = 2
				//					//
				//					// for calculating p(a|bc) in need to know the indices of:
				//					// n(abc) ==> sub2ind(kmer,  0 order  )  >>>> kmer="abc" order=2
				//					// p(c|b) ==> sub2ind(kmer+1,0,order-1)  >>>> kmer=" bc" order=1
				//					// n(ab)  ==> sub2ind(kmer,  0,order-1)  >>>> kmer="abc" order=1
				//					//
				//					// specialities for positions close to the beginning of a sequence (j = 0,1,... j <order)
				//					// let's say my order = 2 and kmer="abc"...
				//					// modelPos = position in the model to start with 0 instead of 1.
				//					// for modelPos=0 --> it is only possible to calculate 0th order models --> order=0
				//					//			      also the kmer needs to be truncated from the front, so a kmer of "abc" will turn to "  c" with kmer+2
				//					// for modelPos=1 --> it is only possible to calculate 1st order models --> order=1
				//					//                    also the kmer needs to be trunctated from the front, so a kmer of "abc" will turn to " bc" with kmer+1
				//					// in a more general way this would look like this:
				//					// n(abc) will fall back to n(bc) for j=1 and n(c) for j=0       ==> sub2ind(kmer, 0, order)
				//					// p(c|b) will stay p(c|b) for j=1 and fall back to xx for j=0  ==> sub2ind(kmer+(order-modelPos), 0, modelPos)
				//					// n(ab)  will fall back to n(b) for j=1 and to xx for j=0      ==> sub2ind(kmer+(order-modelPos), 0, modelPos) //Global::posBg[k];
				//
			}// end learn hyperparameter

		} while( iterate && ( Global::maxEMIterations > round ) );

		if( Global::verbose ){
			printInterpolatedMarkovModel( *motif );
			printf( "\n(after %d iterations)\n", round );
		}

		//		motif->setFactor( motif->getAlpha() );

		if( Global::lastCondsPseudoCounts ){
			motif->resetLastConds();
		}

		for( n=1; n <= N; n++ ){
			free( bg[n] );
			free( likelihoods[n] );
			free( placements[n] );
			free( posterior[n] );
			free( shifts[n] );
		}
		free( bg );
		free( likelihoods );
		free( placements );
		free( posterior );
		free( shifts );

		for( l=0; l < W; l++ ){
			free( lprobs[l] );
		}
		free( lprobs );

		modelCounter++;
	}

	if( Global::saveModels ){
		saveInterpolatedMarkovModels( _models, Global::shortFileName );
	}

	str.str( "" );
	str << Global::outputDirectory << "/" << baseFileName( sequences->name ) << ".logLikelihood";

	FILE *f_likelihood = fopen( str.str().c_str(), "w" );
	fprintf( f_likelihood, "%.8e\n", likelihood );
	fclose( f_likelihood );

	if( Global::testPosSequences ){
		// calculate log odds on training sequences
		score( Global::posSet, Global::logProbs );
	}
	if( Global::testNegSequences ){
		// calculate log odds on background sequences
		score( Global::negSet, Global::logProbs );
	}
}

void EM::score( ss_type sequences, bool logProbs ){

	int N = sequences->nent; // #entities

	if( Global::verbose )
		printf( " _________\n"
				"|         |\n"
				"| SCORING |\n"
				"|_________|\n"
				"\n"
				"on %d sequences (%s)\n", N, sequences->name );

	int b, i, k, l, n;
	int L, LW1, W;

	double **conds;
	double **probs;

	// model placements without N character
	bool** placements;
	// background model shifts
	// due to N character before model placements
	int** shifts;
	int shift;

	double likelihood;
	double likelihoodBg;

	int order;
	int Worder;

	unsigned char* s;
	std::stringstream str;
	FILE* logOddsFile = NULL;

	int nr = 1;
	int modelCounter = 1;

	Motif* motif;
	list<Motif*>::const_iterator iter;
	for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end();
			iter++ ){

		if( Global::nrModels.empty() ){
			nr = modelCounter;
		} else{
			nr = modelCounter;
		}

		if( Global::verbose ){
			printf( " _______\n"
					"|       |\n"
					"| MODEL |\n"
					"|_______|\n"
					"\n"
					"no. %d\n", nr );
		}

		str.str( "" );
		str << Global::outputDirectory << "/" << baseFileName( sequences->name )
		    				<< "-" << nr << ".log" << ( ( logProbs ) ? "Probs" : "Odds" );
		logOddsFile = fopen( str.str().c_str(), "w" );

		motif = *iter;

		motif->logProbs();
		hoNullModel::logProbs();

		b = motif->getFirstMotifColumn();
		W = motif->getMotifLength();

		conds = motif->getConds();
		probs = motif->getPWM();

		order = motif->getOrder();
		Worder = W-order;

		placements = ( bool** )malloc( (N+1)*sizeof( bool* ) );
		shifts = ( int** )malloc( (N+1)*sizeof( int* ) );

		for( n=1; n <= N; n++ ){

			s = sequences->entity[n]->S[0];	// ignore remaining sequences in S
			L = sequences->entity[n]->n;

			placements[n] = ( bool* )calloc( L-W+2, sizeof( bool ) );
			shifts[n] = ( int* )calloc( L-W+2, sizeof( int ) );

			LW1 = L-W+1;
			// precalculate model placements
			// across sequences
			// across sequence positions
			// for subsequences
			for( k=1; k <= LW1; k++ ){
				for( l=0; l < W; l++ ){
					if( !( s[k+l] ) ){
						break; // N character
					}
				}
				if( l == W ){
					// placement possible
					placements[n][k] = true;
					// calculate background model shifts
					if( Global::modelOrderBg > 0 ){
						shift = Global::modelOrderBg;
						l = k-1;
						while( ( l > 0 ) && ( ( k-l ) <= Global::modelOrderBg )
						){
							if( !( s[l] ) ){
								break; // N character
							}
							shift--;
							l--; // back in sequence
						}
						shifts[n][k] = shift;
					}
				}
			}
		}

		// calculate log odds across sequences
		for( n=1; n <= N; n++ ){

			s = sequences->entity[n]->S[0];	// ignore remaining sequences in S
			L = sequences->entity[n]->n;

			if( W > L ){
				fprintf( stderr, "Sequence no. %d shorter than model no. %d\n",
						n, nr );
				fprintf( logOddsFile, "NA\n" );
				continue;
			}

			LW1 = L-W+1;
			// calculate log odds across positions
			for( k=1; k <= LW1; k++ ){

				if( placements[n][k] ){

					// calculate background model probabilities
					likelihoodBg = 0.0;
					// across model positions
					//					for( i=k; (i-k) < W; i++ ){
					//						// using conditional probabilities
					//						likelihoodBg += hoNullModel::getConds()[
					//						                hoNullModel::sub2ind( s,
					//										max( 1, i-Global::modelOrderBg ),
					//										min( i-1, Global::modelOrderBg ) )];
					//					}
					if( !( logProbs ) ){
						shift = shifts[n][k];
						for( i=k; (i-k) < W; i++ ){
							// using conditional probabilities
							likelihoodBg += hoNullModel::getConds()[
							                                        hoNullModel::sub2ind( s,
							                                        		i-Global::modelOrderBg+shift,
							                                        		Global::modelOrderBg-shift )];
							shift = max( 0, shift-1 );
						}
					}

					// calculate model probabilities

					// using max. order probabilities
					likelihood = probs[b+order][motif->sub2ind( s, k, order )];
					// across remaining sequence
					for( i=k+1; (i-k) < Worder; i++ ){
						// using conditional probabilities
						likelihood += conds[b+(i-k)+order][motif->sub2ind(
								s, i, order )];
					}

					// save log odds
					//fprintf( logOddsFile, "%.8e ", likelihood - likelihoodBg
					fprintf( logOddsFile, "%.4e ", likelihood - likelihoodBg
					);
				} else{
					// save placeholder
					fprintf( logOddsFile, "NA " );
				}
			}
			fprintf( logOddsFile, "\n" );
		}
		fclose( logOddsFile );

		motif->expProbs();
		hoNullModel::expProbs();

		for( n=1; n <= N; n++ ){
			free( placements[n] );
		}
		free( placements );

		modelCounter++;
	}
}

void EM::evaluatePWMs( ss_type sequences, bool logProbs ){

	int N = sequences->nent; // #entities

	if( Global::verbose )
		printf( " _________________\n"
				"|                 |\n"
				"| EVALUATING PWMS |\n"
				"|_________________|\n"
				"\n"
				"on %d sequences (%s)\n", N, sequences->name );

	int b, i, k, l, n;
	int L, LW1, W;

	double **probs;

	// model placements without N character
	bool** placements;
	// background model shifts
	// due to N character before model placements
	int** shifts;
	int shift;

	double likelihood;
	double likelihoodBg;

	unsigned char* s;
	std::stringstream str;
	FILE* logOddsFile = NULL;

	int nr = 1;
	int modelCounter = 1;
	list<Motif*>::const_iterator iter;
	for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end();
			iter++ ){

		if( Global::nrModels.empty() ){
			nr = modelCounter;
		} else{
			nr = modelCounter;
		}

		if( Global::verbose ){
			printf( " ___________\n"
					"|           |\n"
					"| PWM MODEL |\n"
					"|___________|\n"
					"\n"
					"no. %d\n", nr );
		}

		str.str( "" );
		str << Global::outputDirectory << "/" << baseFileName( sequences->name )
		    				<< "-XX-" << nr << ".log" << ( ( logProbs ) ? "Probs" : "Odds" );
		logOddsFile = fopen( str.str().c_str(), "w" );

		NullModel::logProbs();

		b = ( *iter )->getFirstMotifColumn();
		W = ( *iter )->getMotifLength();

		probs = ( *iter )->getPWM();

		placements = ( bool** )malloc( (N+1)*sizeof( bool* ) );
		shifts = ( int** )malloc( (N+1)*sizeof( int* ) );

		for( n=1; n <= N; n++ ){

			s = sequences->entity[n]->S[0];	// ignore remaining sequences in S
			L = sequences->entity[n]->n;

			placements[n] = ( bool* )calloc( L-W+2, sizeof( bool ) );
			shifts[n] = ( int* )calloc( L-W+2, sizeof( int ) );

			LW1 = L-W+1;
			// precalculate model placements
			// across sequences
			// across sequence positions
			// for subsequences
			for( k=1; k <= LW1; k++ ){
				for( l=0; l < W; l++ ){
					if( !( s[k+l] ) ){
						break; // N character
					}
				}
				if( l == W ){
					// placement possible
					placements[n][k] = true;
					// calculate background model shifts
					if( Global::order > 0 ){
						shift = Global::order;
						l = k-1;
						while( ( l > 0 ) && ( ( k-l ) <= Global::order )
						){
							if( !( s[l] ) ){
								break; // N character
							}
							shift--;
							l--; // back in sequence
						}
						shifts[n][k] = shift;
					}
				}
			}
		}

		// calculate log odds across sequences
		for( n=1; n <= N; n++ ){

			s = sequences->entity[n]->S[0];	// ignore remaining sequences in S
			L = sequences->entity[n]->n;

			if( W > L ){
				fprintf( stderr, "Sequence no. %d shorter than PWM model no. %d"
						"\n", n, nr );
				fprintf( logOddsFile, "NA\n" );
				continue;
			}

			LW1 = L-W+1;
			// calculate log odds across positions
			for( k=1; k <= LW1; k++ ){

				if( placements[n][k] ){

					// calculate background model probabilities
					likelihoodBg = 0.0;
					// across model positions
					//					for( i=k; (i-k) < W; i++ ){
					//						// using conditional probabilities
					//						likelihoodBg += NullModel::getConds()[
					//						                NullModel::sub2ind( s,
					//										max( 1, i-Global::order ),
					//										min( i-1, Global::order ) )];
					//					}
					if( !( logProbs ) ){
						shift = shifts[n][k];
						for( i=k; (i-k) < W; i++ ){
							// using conditional probabilities
							likelihoodBg += NullModel::getConds()[
							                                      NullModel::sub2ind( s,
							                                    		  i-Global::order+shift,
							                                    		  Global::order-shift )];
							shift = max( 0, shift-1 );
						}
					}

					// calculate model probabilities
					likelihood = 0.0;
					for( l=0; l < W; l++ ){
						likelihood += probs[b+l][s[k+l]];
					}

					// save log odds
					fprintf( logOddsFile, "%.8e ", likelihood - likelihoodBg
					);
				} else{
					// save placeholder
					fprintf( logOddsFile, "NA " );
				}
			}
			fprintf( logOddsFile, "\n" );
		}
		fclose( logOddsFile );

		NullModel::expProbs();

		for( n=1; n <= N; n++ ){
			free( placements[n] );
		}
		free( placements );

		modelCounter++;
	}
}
