#include "backgroundDistribution.h"
#include "Globals.h"
#include "LogTable.h"
#include "StartModels.h"
#include "utils.h"
#include "em/EM.h"
#include "em/hoNullModel.h"
#include "em/hoUtils.h"
#include "nucleotides/conservationScores.h"
#include "refinementPhase/MotifContainer.h"
#include "refinementPhase/Refine_Motifs.h"

int main(int argc, char *argv[])
{
    long time1 = time( NULL );

	fprintf( stderr, "\n" );
	fprintf( stderr, "================\n" );
	fprintf( stderr, "=   XX MOTIF   =\n" );
	fprintf( stderr, "================\n" );

    // initialize global variables
	Global G( argc, argv );

	// pre-calculate logarithms
	LogTable L;

	// positive sequence set statistics
	printf( "posSet statistics: %d, maxMultSeq: %d, minLength: %d, maxLength: "
			"%d, total: %d\n", Global::posSet->nent,
			Global::posSet->max_MultSeq, Global::posSet->min_leng,
			Global::posSet->max_leng,
			Global::posSet->total[Global::posSet->nent]);
	// negative sequence set statistics
	if( Global::negSet != NULL ){
		printf("negSet statistics: %d, maxMultSeq: %d, minLength: %d, "
				"maxLength: %d, total: %d\n", Global::negSet->nent,
				Global::negSet->max_MultSeq, Global::negSet->min_leng,
				Global::negSet->max_leng,
				Global::negSet->total[Global::negSet->nent]);
	}


	MotifContainer startModels;
	StartModels seeds;

	if( Global::em ){

		// initialize model probabilities
		if( Global::bindingSiteFile == NULL && Global::markovModelFile == NULL
			){
			// initialize from XXmotif results
			if( Global::profFile == NULL && Global::startMotif == NULL ){
				seeds.findInitialMotifs( startModels, 1e-3 );
			} else{
				seeds.initStartMotif( startModels );
			}
			MotifRefinement refinement( startModels );
			refinement.start();
		} else{
			// initialize from binding sites or Markov model file
			seeds.initStartMotif( startModels );
		}

		if( !( Global::noExpectationMaximizationPhase ) ){
			if( startModels.getMotifs().size() > 0 ){
				// EM phase
				EM em( startModels );
				em.go();
				if( !( Global::testSet.empty() ) ){
					// calculate log odds on test sequences
					std::list<ss_type>::const_iterator iter;
					for( iter=Global::testSet.begin(); iter != Global::testSet.end(); iter++ ){
						em.score( *iter, Global::logProbs );
					}
				}
			} else{
				fprintf( stderr, "No higher-order model to process.\n" );
			}
		} else{
			if( startModels.getMotifs().size() > 0 ){
				if( !( Global::testSet.empty() ) ){
					// no EM phase
					EM em( startModels );
					// calculate log odds on test sequences
					std::list<ss_type>::const_iterator iter;
					for( iter=Global::testSet.begin(); iter != Global::testSet.end(); iter++ ){
						em.score( *iter, Global::logProbs );
					}
				}
			}
		}
	} else{ // XXmotif or MadonaPro


		if( Global::profFile == NULL && Global::startMotif == NULL ){
			seeds.findInitialMotifs( startModels, 1e-3 );
		} else{
			seeds.initStartMotif( startModels );
		}

		MotifRefinement refinement( startModels );
		refinement.start();
	}
	printf( "\n" );

	fprintf( stderr, "\n------------ time: %ld seconds (%0.2f minutes) "
			"------------\n", time( NULL ) - time1, ( float )( time( NULL ) -
			time1 )/60.0 );

	return 0;	
}
