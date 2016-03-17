#include "MotifContainer.h"
#include "../ProgressBar.h"

MotifContainer::~MotifContainer(){
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end();){
		delete *it;
		it = _startModels.erase(it);
	}
}

void MotifContainer::clear(){
	cerr << "clear MotifContainer" << endl;
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end();){
		delete *it;
		it = _startModels.erase(it);
	}
}

void MotifContainer::update_Pval_and_sort(double pseudo){
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end();){
		(*it)->updatePval((*it)->isTracked());

		if((*it)->isTracked()) cerr << "\n+++++++++++ Tracked Motif: " << **it;

		if((*it)->getStartPosList().size() == 0){
			if((*it)->isTracked()) cerr << "\n+++++++++++ Tracked Motif has no seeds => filtered out " << **it;
			delete *it;
			it = _startModels.erase(it);
			continue;
		}
		it++;
	}
	_startModels.sort(Motif::cmpPtr);
}

void MotifContainer::updatePWM(double pseudo){
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end(); it++){
		(*it)->updatePWM_OOL_initial(pseudo);
	}
}

void MotifContainer::setBindingSites(const ThresholdChecker &pValThreshold){
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end(); it++){
		(*it)->setBindingSites(pValThreshold);
	}
}

void MotifContainer::filter(const int _minimumMotifs, const int _maximumMotifs, const double _pValThreshold){
	int motNb = 0;
	int filteredMotifs = 0;
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end(); ){
		if((motNb < _minimumMotifs || (*it)->getPval() < _pValThreshold) && motNb < _maximumMotifs){
			motNb++;
			it++;
		}else{
			motNb++;
			if((*it)->isTracked()){
				fprintf(stderr, "+++++++ Tracked Motif is filtered out\n" \
						"+++	Thresholds: minimumMotfs: %d, E-Value: %f\n" \
						"+++	Motif: motNb: %d, E-Value: %f\n\n", _minimumMotifs, _pValThreshold, motNb, (*it)->getPval() );
			}
			delete *it;
			it = _startModels.erase(it);
			filteredMotifs++;
		}
	}
	//if(filteredMotifs != 0){
	//	fprintf(stderr, "\n%d motifs above E-value threshold of %.2e => filtered out\n", filteredMotifs, exp(_pValThreshold));
	//}
}

void MotifContainer::removeHomologousMotifs() {
	for (list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end(); ) {
		if (static_cast<int>((*it)->getMotifColumns().size()) == Global::maxMatchPositions) {
			bool all_equal = true;
			std::vector<int> residues;
			const StartPosContainer& sp = (*it)->getStartPosList();
			const motif_columns_type& cols = (*it)->getMotifColumns();
			for (StartPosContainer::const_iterator pit = sp.begin(); pit != sp.end(); ++pit) {
				int num = 0;
				if (residues.empty()) {
					for (motif_columns_type::const_iterator cit = cols.begin(); cit != cols.end(); ++cit, ++num) {
						const int aa = Global::posSet->entity[pit->seq]->S[0][pit->pos+*cit-1];
						residues.push_back(aa);
					}
				} else {
					for (motif_columns_type::const_iterator cit = cols.begin(); cit != cols.end(); ++cit, ++num) {
						const int aa = Global::posSet->entity[pit->seq]->S[0][pit->pos+*cit-1];
						if (aa != residues[num]) {
							all_equal = false;
							break;
						}
					}
				}
			}
			if (all_equal) {
				cerr << "WARNING: Conservation filter: Removing " << (**it).getIUPACString() << endl;
				it = _startModels.erase(it);
			} else {
				++it;
			}
		} else {
			++it;
		}
	}
}

void MotifContainer::removeRedundantMotifs(){
	list<Motif*>::iterator next = _startModels.begin();
	if(next == _startModels.end()) return;
	for(list<Motif*>::iterator it = _startModels.begin(); ++next != _startModels.end(); ){

		int identical = 0;
		/* check whether two motifs are completely identical */
		if((*it)->getStartPosList().size() == (*next)->getStartPosList().size()){
			identical = 0;

			StartPosContainer& sp1 = (*it)->getStartPosList();
			StartPosContainer& sp2 = (*next)->getStartPosList();

			StartPosContainer::iterator it1 = sp1.begin();
			StartPosContainer::iterator it2 = sp2.begin();

			while(*it1 == *it2){
				it1++;
				it2++;
				if(it1 == sp1.end() && it2 == sp2.end()){
					identical = 1;
					break;
				}
			}
		}
		if(identical){
			if((*it)->isTracked()){	(*next)->setTracked();	}
			delete *it;
			it = _startModels.erase(it);
		}else{
			it++;
		}
	}
}

void MotifContainer::merge(bool iterativePhase){
	/* merge PWMs with each other */
	int round = 1;
	bool loop = false;
	bool removeOverlappingMotifs = false;
	std::cout << " Merge Motifs: " << std::flush;
	do{
		loop = false;
		/* set iterator to first mergable motif ( second one ) */
		int count = 0;
		list<Motif*>::iterator it = _startModels.begin(); it++; count++;

		/* test all motifs whether they can be merged to the already merged motifs */
		for(int i=0; it != _startModels.end(); i++){
		    if(i%500 == 0) std::cout << "." << std::flush;
		    int del = 0;
		    int count_merged = 0;
			for(list<Motif*>::iterator merged_it = _startModels.begin(); merged_it != it; merged_it++, count_merged++){
				/* test whether the already accepted Motif can be merged to the newly tested motif */

				int merge = -10000;
				int mergeRevcomp = 0;
				/* merge before iterative phase only if two motifs have the same type */

				//if(iterativePhase || abs((*merged_it)->getMotifLength() - (*it)->getMotifLength()) < 3){
					merge = MotifComparison::calcBestMergingDist(*merged_it, *it, mergeRevcomp);
				//}

				/* merge with given merging distance*/
				if(merge != -10000){
					if(MotifComparison::DEBUG){
						fprintf(stderr, "\nmerge model %d with %d : %d\n", count_merged, count, merge);
						cerr << **merged_it << **it;

					}
					del = 1;
					//fprintf(stderr, "start merge\n");
					//cerr << **merged_it << **it;
					MotifComparison::mergePWMwithPWM(*merged_it, *it, merge, mergeRevcomp);
					//fprintf(stderr, "end merge\n");
					(*merged_it)->resetConverged();
					/* merge only with one other motif */
				}else{
					if(removeOverlappingMotifs){
						if(MotifComparison::DEBUG)	cerr <<"remove overlapping motifs: model " << count_merged << " and " << count << " : " << merge << endl << **merged_it << **it;
						MotifComparison::removeOverlappingMotifOccurrencies(*merged_it, *it, (-1)*merge);
						(*merged_it)->resetConverged();
						(*it)->resetConverged();
						if(MotifComparison::DEBUG){	cerr << "after removal: " << endl << **merged_it << **it;	exit(-1); }
					}
				}
				//exit(-1);
			}
			if(del){
				delete *it;
				it = _startModels.erase(it);
				/* change occurred: do another round of merging */
				loop = true;
			}else{
				it++;
			}
			count++;
		}

		_startModels.sort(Motif::cmpPtr);

		cout << "\tround " << round++ << ": " << _startModels.size() << "\t";
	}while(loop);
	cout << endl;

	updatePWM(Global::pseudo);
}

void MotifContainer::iterate(int maximizeMotifLength, bool lastIteration){
	IterateMotif iterateMotif;

	int motifNb = 1;
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end(); it++, motifNb++){
		if(!Global::batch_mode){
			cout << "\riterate Motifs: " << progressBar(motifNb, static_cast<int>(_startModels.size()), 40) << std::flush;
		}

		Motif* m = *it;

		double oldPval=1000, bestPval=1000;
		int loop = 1;

		if(Global::DEBUG || m->isTracked()) cerr << endl << "start Motif" << endl << *m;
		bool keepPosSetSize = true;
		do{
			if(m->isConverged()) break;
			oldPval = bestPval;
			iterateMotif.optimizeMotif(m, keepPosSetSize, false);
			if (Global::usePositionalProbs) {
				/* update enriched region:
				   use best matrix matches without positional prior and calculate new enriched region
				   with all matches within a region which is twice as large as the old one
				   => if region stays small, there is really an enrichment */
				m->updateEnrichedRegion(StartPosUpdater::getInstance().getBestStartPos());
			}

			if(Global::DEBUG || m->isTracked()) cerr << endl << "optimize plus" << endl << *m;

			/* find the beset motif length if the length has not already converged */
			bool lengthCorrection = maximizeMotifLength; // && !length_converged;

			iterateMotif.maximizeMotifLength(m, lengthCorrection, false);

			bestPval = m->getPval();
			if(Global::DEBUG || m->isTracked()) cerr << "after length optimization:" << *m;

			/* check whether the motif has already converged */
			if(loop > 0 && bestPval >= oldPval){
				loop=0;
			}

		}while(loop > 0);

		/* sort start Positions */
		m->getStartPosList().sort(StartPos::cmpStartPos);

		m->setConverged();

		/* print debug information */
		if(Global::DEBUG) cout << *m;
		if(m->isTracked()) cerr << endl << "+++++++++++ Tracked Motif:" << *m;
	}
}

void MotifContainer::sort_and_filter(double filterThreshold){
	_startModels.sort(Motif::cmpPtr);

	/* only accept motifs with a pValue better then _filterThreshold for the next round */
	int count = 0;
	for(list<Motif*>::iterator it = _startModels.begin(); it != _startModels.end();){
		if((count >= 20 && (*it)->getPval() > filterThreshold)){
			if((*it)->isTracked()){
				cerr << endl <<
						"+++++++++++ Tracked Motif is filtered out:" << endl;
				cerr << "++\t\tmotifNb: " << count+1 << "\tmin Threshold: 20" << endl;
				cerr << "++\t\tpVal: " << (*it)->getPval() << "\tmax Threshold: " << filterThreshold << endl;
				cerr << "++\t\tsites: " << (*it)->getTotalSites() << "\tmin Threshold: " << 1 << endl << endl;
			}
			delete *it;
			it = _startModels.erase(it);
		}else{
			count++;
			it++;
		}
	}
}

void MotifContainer::sortByPosPval(){
	_startModels.sort(Motif::cmpPosPtr);
}

void MotifContainer::sort() {
	_startModels.sort(Motif::cmpPtr);
}

/*
 * filter by model number in ranking
 */
void MotifContainer::filter( std::vector<int> nrModels ){

	if( nrModels.size() > 0 ){

		int i = 0;
		int modelCounter = 1;
		int filterCounter = 0;

		for( list<Motif*>::iterator iter = _startModels.begin(); iter != _startModels.end(); ){
			if( i < static_cast<int>( nrModels.size() ) ){
				if( modelCounter == nrModels.at( i ) ){
					i++;
					iter++;
				} else{
					if( ( *iter )->isTracked() ){
						fprintf( stderr, "Tracked motif fails to pass filter\n" );
					}
					delete *iter;
					iter = _startModels.erase( iter );
					filterCounter++;
				}

			} else{
				if( ( *iter )->isTracked() ){
					fprintf( stderr, "Tracked motif fails to pass filter\n" );
				}
				delete *iter;
				iter = _startModels.erase( iter );
				filterCounter++;
			}
			modelCounter++;
		}
	}
}

/*
 * filter min. percentage of sequences containing a binding site instance
 * in addition to
 * min. number of models
 * max. number of models
 * max. p-value of models
 * */
void MotifContainer::filter( const int minModels, const int maxModels,
		                     const double maxPvalue, const float minOccurrence ){

	int lastSeq;
	int uniqueSeqs;
	float occurrence;

	int modelCounter = 0;
	int filterCounter = 0;

	for( list<Motif*>::iterator iter = _startModels.begin();
		 iter != _startModels.end(); ){

		if( Global::multipleOccurrence ){

			lastSeq = 0;
			uniqueSeqs = 0;
			StartPosContainer sites = ( *iter )->getStartPosList();
			for( StartPosContainer::iterator siteiter=sites.begin(); siteiter !=
				 sites.end(); siteiter++ ){ // binding site iterator

				if( siteiter->seq != lastSeq ){
					++uniqueSeqs;
				}
				lastSeq = siteiter->seq;
			}
		} else{

			uniqueSeqs = ( *iter )->getTotalSites();
		}
		occurrence = static_cast<float>( uniqueSeqs ) /
				     static_cast<float>( ( *iter )->getPosSetSize() );

		if( ( modelCounter < minModels ||
			( (*iter)->getPval() <= maxPvalue && occurrence >= minOccurrence ) )
			&& modelCounter < maxModels ){

			modelCounter++;
			iter++;
		}else{
			if( ( *iter )->isTracked() ){
				fprintf( stderr, "Tracked motif fails to pass filter\n" );
			}
			delete *iter;
			iter = _startModels.erase( iter );
			filterCounter++;
		}
	}
}

/*
 * initialize higher-order models from XXmotif results
 * */
void MotifContainer::initHoMotifsWithStartPos(){

	for( list<Motif*>::iterator iter = _startModels.begin(); iter != _startModels.end(); iter++ ){
		( *iter )->initHoMotifWithStartPos();
	}
}

/*
 * calculate the percentage of positive sequences containing a binding site
 * instance used as model-specific specificity factor in EM calculations
 */
void MotifContainer::setSpecificityFactors(){

	// last binding site sequence number
	int lastSeq;
	// number of positive sequences with one (or more) binding site instance(s)
	int uniqueSeqs;
	// model-specific specificity factor
	float q;
	// binding site sequences
	StartPosContainer sites;

	for( list<Motif*>::iterator modeliter = _startModels.begin(); modeliter !=
		 _startModels.end(); modeliter++ ){ // model iterator

		// binding site sequences
		sites = ( *modeliter )->getStartPosList();

		if( Global::multipleOccurrence ){ // one or more instances per sequence
			
			lastSeq = 0;
			uniqueSeqs = 0;
			for( StartPosContainer::iterator siteiter=sites.begin(); siteiter !=
				 sites.end(); siteiter++ ){ // binding site iterator

				if( siteiter->seq != lastSeq ){
					++uniqueSeqs;
				}
				lastSeq = siteiter->seq;
			}
		} else{ // max. one instance per sequence
			uniqueSeqs = ( *modeliter )->getTotalSites();
		}

		// calculate model-specific specificity factor
		q = static_cast<float>( uniqueSeqs ) /
			static_cast<float>( ( *modeliter )->getPosSetSize() );
		// set model-specific specificity factor
		( *modeliter )->setSpecificityFactor( std::min( q, Global::qmax ) );
	}
}
