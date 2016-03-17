#include "Refine_Motifs.h"
#include "../output.h"
//#include "../aminoacids/madonaPro.h"
#include "../AbstractKmer.h"
#include "../Globals.h"
#include "../SmallKmer.h"
#include "../elongationPhase/Kmer.h"
#include "../elongationPhase/Match.h"
#include "../elongationPhase/elongationCandidates.h"
#include "../memoryPool/pool_alloc.h"

#include "../em/EM.h"
#include "../em/hoUtils.h"

/*find best Motifs by iterating till convergence */
void MotifRefinement::start(){
	int oldMergedMotifsNb = _models.getMotifNb();

	printf("Number of motifs: %d\n\n", _models.getMotifNb());

	if(!Global::noRefinementPhase){
		printf("\n\tITERATIVELY REFINE PWMS\n\n");

		bool lastIteration = false;
		int round = 1;
		do{
			if(round != 1 && _models.getMotifNb() == oldMergedMotifsNb) lastIteration = true;
			/* find best motif by iterating them */
			if (Global::maxIterations > 0){
				_models.iterate(Global::maximizeMotifLength, lastIteration);
			}
			if(_models.getMotifNb() == 0 || Global::maxIterations<=1) lastIteration = true;

			/* sort found motifs*/
			_models.sort_and_filter(100);

			if(!lastIteration){
				oldMergedMotifsNb = _models.getMotifNb();
				/* merge motifs and remove overlapping motifs */
				//cerr << endl << "merge models" << endl;
				_models.merge(true);
				//cerr << endl << "merged" << endl;
				//_models.merge(false);
			}
			round++;
		}while(!lastIteration);
	}

	/* Rescale E-values by taking square root.
	 * (Improved correspondance between theoretical and reported
	 * E-values on random sequences)
	 */
	const list<Motif*> &motifs = _models.getMotifs();
	for( list<Motif*>::const_iterator mit = motifs.begin(); mit != motifs.end();
		 ++mit ){
		( *mit )->setPval( 0.5 * ( *mit )->getPval() );
	}

	int maxFilter = std::numeric_limits<int>::max();
	if( Global::noRefinementPhase ){
		maxFilter = 10;
	}

	if( Global::positionalProbsRanking ){
		_models.sortByPosPval();
		_models.filter( Global::minModels, Global::maxModels,
				        Global::pValueThreshold, Global::minOccurrence );
	} else{
		_models.filter( 5, maxFilter, Global::finalFilterThreshold );
	}


	_models.setBindingSites(Global::instanceThreshold);

	if (Global::removeHomology) {
		_models.removeHomologousMotifs();
	}


	/* write results into files */
	std::cout << std::endl << "********** Final results **********" << std::endl << std::endl;

	Output::printOutput(_models, 10);

	if( Global::em ){

		if( !( Global::positionalProbsRanking ) ){
			_models.filter( Global::minModels, Global::maxModels,
					        Global::pValueThreshold, Global::minOccurrence );
		}
		if( !( Global::nrModels.empty() ) ){
			_models.filter( Global::nrModels );
		}
		if( _models.getMotifs().empty() ){
			fprintf( stderr, "No motifs pass filter\n" );
			exit(1);
		}

		if( Global::msq ){
			// calculate model-specific specificity factors
			_models.setSpecificityFactors();
		}

		if( Global::evaluatePWMs ){

			// evaluate PWM model(s)
			EM em( _models );

			if( Global::testPosSequences ){
				em.evaluatePWMs( Global::posSet, Global::logProbs );
			}
			if( Global::testNegSequences ){
				em.evaluatePWMs( Global::negSet, Global::logProbs );
			}
			if( !( Global::testSet.empty() ) ){
				std::list<ss_type>::const_iterator iter;
				for( iter=Global::testSet.begin(); iter !=
					 Global::testSet.end(); iter++ ){
					em.evaluatePWMs( *iter, Global::logProbs );
				}
			}
		}

		// initialize higher-order models from XXmotif results
		_models.initHoMotifsWithStartPos();

		if( Global::saveInitModels ){
			std::stringstream str;
			str.str( "" );
			str << Global::shortFileName << "-init";
			saveInterpolatedMarkovModels( _models, str.str().c_str() );
		}
	} else if( Global::evaluatePWMs ){

		if( !( Global::positionalProbsRanking ) ){
			_models.filter( Global::minModels, Global::maxModels,
					        Global::pValueThreshold, Global::minOccurrence );
		}
//		printf( "Motifs:\n" );
//		for( list<Motif*>::const_iterator iter = _models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){
//			printf( "%s\n", ( *iter )->getIUPACString().c_str() );
//		}
		if( !( Global::nrModels.empty() ) ){
			_models.filter( Global::nrModels );
		}
//		printf( "Motifs after filtering:\n" );
//		for( list<Motif*>::const_iterator iter = _models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){
//			printf( "%s\n", ( *iter )->getIUPACString().c_str() );
//		}
		if( _models.getMotifs().empty() ){
			fprintf( stderr, "No motifs pass filter\n" );
		} else{
			// evaluate PWM model(s)
			EM em( _models );

			if( Global::testPosSequences ){
				em.evaluatePWMs( Global::posSet, Global::logProbs );
			}
			if( Global::testNegSequences ){
				em.evaluatePWMs( Global::negSet, Global::logProbs );
			}
			if( !( Global::testSet.empty() ) ){
				std::list<ss_type>::const_iterator iter;
				for( iter=Global::testSet.begin(); iter !=
					 Global::testSet.end(); iter++ ){
					em.evaluatePWMs( *iter, Global::logProbs );
				}
			}
		}
	}
}
