#include "Motif.h"
#include "Compare_Motifs.h"
#include "Iterate_Motif.h"
#include "StartPosUpdater-inl.h"
#include "../Globals.h"
#include "../pValCalculation.h"
#include "../branch_and_bound-inl.h"

#include "../em/hoUtils.h"
#include "../em/hoNullModel.h"

#include <stdbool.h>
#include <stdio.h>

#include <limits>

Motif::Motif(int length) : _length(length){

	_converged = false;
	_countStartPos.counts = ( int* )calloc( Global::posSet->max_leng+1, sizeof( int ) );
	_countStartPos.size = 0;
	_posSetSize = Global::posSet->nent;
	_tracked = false;

	_order = Global::modelOrder;
	_alpha = Global::alpha;
	_q = Global::q;

	if( Global::em && ( Global::bindingSiteFile != NULL ||
			Global::markovModelFile != NULL ) ){


		_coeffs = ( int* )calloc( _order+2, sizeof( int ) );
		_offsets = ( int* )calloc( _order+2, sizeof( int ) );

		int k;
		for( k=0; k < _order+2; k++ ){
			_coeffs[k] = static_cast<int>( pow( nAlpha( Global::A ), k ) );
			_offsets[k] = k ? _offsets[k-1]+_coeffs[k] : _coeffs[k];
		}
		_fields = _offsets[k-1];

		_counts = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_counts[i] = ( double* )calloc( _fields, sizeof( double ) );
		}
		_countsx = ( double** )malloc( ( _length )*sizeof( double* ) );
		for( int i=0; i < _length; i++ ){
			_countsx[i] = ( double* )calloc( _offsets[_order], sizeof( double )
			);
		}

		_conds = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_conds[i] = ( double* )calloc( _fields, sizeof( double ) );
		}
		_lastConds = NULL;

	} else{

		_coeffs = NULL;
		_offsets = NULL;
		_fields = nAlpha( Global::A ) + 1;

		_counts = NULL;
		_countsx = NULL;

		_conds = NULL;
		_lastConds = NULL;
	}

	_pwm0 = NULL;

	_pwm = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
	for( int i=0; i <= _length; i++ ){
		_pwm[i] = ( double* )calloc( _fields, sizeof( double ) );
	}
}

Motif::Motif(const std::shared_ptr<Kmer>& skr, int length){

	_length = length;
	_converged = false;

	const AbstractKmer* kmer = skr->getKmer();
	_tracked = Global::isTracked( kmer->bestNucString());

	int pos = ( _length - kmer->length() ) / 2 + 1;
	_motifColumns.push_back( pos++ );
	for( int i=0; i < kmer->numMatches()-1; i++ ){
		pos += kmer->gapsAfter( i ); /* omitted (gap) positions */
		_motifColumns.push_back( pos++ );
	} // e.g. 7 | 8 | 9 | 20 | 21 | 22 |

	_enrichment.max = skr->enrichment.max - getFirstMotifColumn() + 1;
	_enrichment.startRegion = skr->enrichment.startRegion - getFirstMotifColumn() + 1;
	_enrichment.endRegion = skr->enrichment.endRegion - getFirstMotifColumn() + 1;

	_pVal = skr->p_set;

	_order = Global::modelOrder;
	_alpha = Global::alpha;
	_q = Global::q;

	if( Global::em && ( Global::bindingSiteFile != NULL ||
			Global::markovModelFile != NULL ) ){


		_coeffs = ( int* )calloc( _order+2, sizeof( int ) );
		_offsets = ( int* )calloc( _order+2, sizeof( int ) );

		int k;
		for( k=0; k < _order+2; k++ ){
			_coeffs[k] = static_cast<int>( pow( nAlpha( Global::A ), k ) );
			_offsets[k] = k ? _offsets[k-1]+_coeffs[k] : _coeffs[k];
		}
		_fields = _offsets[k-1];

		_counts = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_counts[i] = ( double* )calloc( _fields, sizeof( double ) );
		}
		_countsx = ( double** )malloc( ( _length )*sizeof( double* ) );
		for( int i=0; i < _length; i++ ){
			_countsx[i] = ( double* )calloc( _offsets[_order], sizeof( double )
			);
		}

		_conds = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_conds[i] = ( double* )calloc( _fields, sizeof( double ) );
		}
		_lastConds = NULL;

	} else{

		_coeffs = NULL;
		_offsets = NULL;
		_fields = nAlpha( Global::A ) + 1;

		_counts = NULL;
		_countsx = NULL;

		_conds = NULL;
		_lastConds = NULL;
	}

	_pwm0 = NULL;

	_pwm = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
	for( int i=0; i <= _length; i++ ){
		_pwm[i] = ( double* )calloc( _fields, sizeof( double ) );
	}

	fill_startPos_with_IUPAC_matches( skr );

	int first = getFirstMotifColumn();
	int end = getLastMotifColumn();

	_motifColumns.clear();
	for( int i=first; i <= end; i++ )
		_motifColumns.push_back( i );
	//	 // e.g. 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 |

	/* Update motif enrichment */

	_countStartPos.counts = ( int* )calloc( Global::posSet->max_leng+1, sizeof( int ) );
	_countStartPos.size = 0;

	if( Global::usePositionalProbs ){

		const int startRegion = getEnrichment().startRegion + getFirstMotifColumn() - 1;
		const int endRegion = getEnrichment().endRegion + getFirstMotifColumn() - 1;
		const int region_ext = static_cast<int>( round( 0.5*( _enrichment.endRegion - _enrichment.startRegion + 1 ) ) );

		for( StartPosContainer::iterator it = _sites.begin(); it!= _sites.end(); it++ ){
			int startPos = it->pos + getFirstMotifColumn()-1;
			if(Global::revcomp && startPos > Global::posSet->max_leng / 2){
				startPos = Global::posSet->max_leng - (startPos + getMotifLength() - 1);
			}
			if( startPos >= startRegion - region_ext && startPos <= endRegion + region_ext ){
				_countStartPos.counts[startPos]++;
				_countStartPos.size++;
			}
		}
		int maxRegion = Global::posSet->max_leng - getMotifLength() + 1;
		if(Global::revcomp){
			maxRegion = (Global::posSet->max_leng - 1) / 2 - getMotifLength() + 1;
		}
		MotifRegion motifRegion(1, maxRegion);
		region r = motifRegion.getRegion(_countStartPos.counts, _countStartPos.size);
		r.max-=getFirstMotifColumn()-1; r.startRegion-=getFirstMotifColumn()-1; r.endRegion-= getFirstMotifColumn()-1;

		setEnrichment(r);
	} else{
		region r;
		r.set = 0;
		r.max = 0;
		r.startRegion = -getFirstMotifColumn();
		r.endRegion = Global::posSet->max_leng - getMotifLength() + getFirstMotifColumn();
		setEnrichment( r );
	}

	/*
	if( isTracked() ){
		updatePWM_OOL_initial( Global::pseudo );
		cerr << *this;
		updatePval( false);
		cerr << *this;
	}
	 */
}

Motif::Motif(const Motif& other){
	_length = other._length;
	_enrichment = other._enrichment;
	_motifColumns = other._motifColumns;
	_pVal = other._pVal;
	_converged = false;
	_tracked = other._tracked;

	_countStartPos.counts = (int*)calloc(Global::posSet->max_leng+1, sizeof(int));
	for(int i=0; i <= Global::posSet->max_leng; i++){
		_countStartPos.counts[i] = other._countStartPos.counts[i];
	}
	_countStartPos.size = other._countStartPos.size;
	_sites = other._sites;

	_posSetSize = other._posSetSize;

	_order = other._order;

	_alpha = other._alpha;

	_q = other._q;

	// this was added by me to overcome warnings on uninitialized parameter....
	_worstScore = other._worstScore;
	_posPVal = other._posPVal;
	// end

	// used during alpha optimization
	_Qfunc_p1 = 0.0;
	_Qfunc_p2 = 0.0;
	_Qfunc_p3 = 0.0;
	_Qfunc_p4 = 0.0;
    _Qfunc_grad = 0.0;
	_Grad_p1 = 0.0;
	_Grad_p2 = 0.0;
	_Grad_p3 = 0.0;

	_zeroOrder = 0.0;
	_pseudoOrder = 0.0;
	// end

	if( other._coeffs != NULL ){

		_coeffs = ( int* )calloc( _order+2, sizeof( int ) );

		int k;
		for( k=0; k < _order+2; k++ ){
			_coeffs[k] = other._coeffs[k];
		}
	} else{

		_coeffs = NULL;
	}

	if( other._offsets != NULL ){

		_offsets = ( int* )calloc( _order+2, sizeof( int ) );

		int k;
		for( k=0; k < _order+2; k++ ){
			_offsets[k] = other._offsets[k];
		}
	} else{

		_offsets = NULL;
	}
	_fields = other._fields;

	if( other._counts != NULL ){

		_counts = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_counts[i] = ( double* )calloc( _fields, sizeof( double ) );
			for( int k=0; k < _fields; k++ ){
				_counts[i][k] = other._counts[i][k];
			}
		}
	} else{

		_counts = NULL;
	}

	if( other._countsx != NULL ){

		_countsx = ( double** )malloc( ( _length )*sizeof( double* ) );
		for( int i=0; i < _length; i++ ){
			_countsx[i] = ( double* )calloc( _offsets[_order], sizeof( double )
			);
			for( int k=0; k < _offsets[_order]; k++ ){
				_countsx[i][k] = other._countsx[i][k];
			}
		}
	} else{

		_countsx = NULL;
	}

	if( other._conds != NULL ){

		_conds = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_conds[i] = ( double* )calloc( _fields, sizeof( double ) );
			for( int k=0; k < _fields; k++ ){
				_conds[i][k] = other._conds[i][k];
			}
		}
	} else{

		_conds = NULL;
	}

	if( other._lastConds != NULL ){

		_lastConds = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_lastConds[i] = ( double* )calloc( _fields, sizeof( double ) );
			for( int k=0; k < _fields; k++ ){
				_lastConds[i][k] = other._lastConds[i][k];
			}
		}
	} else{

		_lastConds = NULL;
	}

	if( other._pwm0 != NULL ){

		_pwm0 = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
		for( int i=0; i <= _length; i++ ){
			_pwm0[i] = ( double* )calloc( nAlpha( Global::A )+1, sizeof( double ) );
			for( int k=0; k < ( nAlpha( Global::A )+1 ); k++ ){
				_pwm0[i][k] = other._pwm0[i][k];
			}
		}
	} else{

		_pwm0 = NULL;
	}

	_pwm = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
	for( int i=0; i <= _length; i++ ){
		_pwm[i] = ( double* )calloc( _fields, sizeof( double ) );
		for( int j=0; j<_fields; j++ ){
			_pwm[i][j] = other._pwm[i][j];
		}
	}
}

Motif::~Motif(){

	for( int i=0; i <= _length; i++ ){
		free( _pwm[i] );
	}
	free( _pwm );
	free( _countStartPos.counts );

	if( Global::em ){

		if( _counts != NULL ){
			for( int i=0; i <= _length; i++ ){
				free( _counts[i] );
			}
		}
		free( _counts );

		if( _countsx != NULL ){
			for( int i=0; i < _length; i++ ){
				free( _countsx[i] );
			}
		}
		free( _countsx );

		if( _conds != NULL ){
			for( int i=0; i <= _length; i++ ){
				free( _conds[i] );
			}
		}
		free( _conds );

		if( _lastConds != NULL ){
			for( int i=0; i <= _length; i++ ){
				free( _lastConds[i] );
			}
		}
		free( _lastConds );

		if( _coeffs != NULL ){
			free( _coeffs );
		}

		if( _offsets != NULL ){
			free( _offsets );
		}

		if( _pwm0 != NULL ){
			for( int i=0; i <= _length; i++ ){
				free( _pwm0[i] );
			}
		}
		free( _pwm0 );
	}
}

void Motif::updatePWM(double pseudo, bool maximizeMotifLength){
	updateCountsPWM(maximizeMotifLength);
	normalizeLOG(pseudo, maximizeMotifLength);
}

void Motif::filter_region(){
	for (StartPosContainer::iterator it=_sites.begin(); it!=_sites.end();){
		if(it->pos < _enrichment.startRegion || it->pos > _enrichment.endRegion){
			it = _sites.erase(it);
		}else{ it++; }
	}
}

void Motif::filter_fitToSequence(){
	for (StartPosContainer::iterator it=_sites.begin(); it!=_sites.end();){
		if(it->pos + getFirstMotifColumn() - 1 < 1 || it->pos + getLastMotifColumn() - 1 > Global::posSet->entity[it->seq]->n){
			it = _sites.erase(it);
		}else if(it->seq > getPosSetSize()){
			it = _sites.erase(it);
		}
		else{ it++; }
	}
}

void Motif::filter_oops_updatePWM(){
	/* update counts for pwms */
	updateCountsPWM();

	StartPosContainer::iterator it_next = _sites.begin(); it_next++;
	for (StartPosContainer::iterator it = _sites.begin(); it_next != _sites.end(); it++, it_next++){

		/* check whether there are is one sequence with two motifs */
		while(it_next != _sites.end() && it->seq == it_next->seq){
			int seqLen = Global::posSet->entity[it->seq]->n;
			double bestScore = 0, score = 0;

			/* decide which binding site should be removed */
			uint8_t* S = Global::posSet->entity[it->seq]->S[0];

			for(motif_columns_type::const_iterator it_column = _motifColumns.begin(); it_column != _motifColumns.end(); it_column++){
				int pos = it->pos + *it_column - 1;
				int pos_next = it_next->pos + *it_column - 1;
				int base, base2;
				if(pos > 0 && pos <= seqLen)
					base = S[pos];
				else base = 0;

				if(pos_next > 0 && pos_next <= seqLen)
					base2 = S[pos_next];
				else base2 = 0;

				if(base == 0){	bestScore += LogTable::LOG1_1000;	continue;	}
				if(base2== 0){	score += LogTable::LOG1_1000;		continue;	}

				bestScore += fast_log(static_cast<float>(_pwm[*it_column][base] / _pwm[*it_column][0]));
				score += fast_log(static_cast<float>(_pwm[*it_column][base2] / _pwm[*it_column][0]));
			}

			/* a better motif in the current sequence has been found => remove current motif */
			if(score > bestScore){
				bestScore = score;
				int pos = it->pos+getFirstMotifColumn()-1;
				for(int k=getFirstMotifColumn(); k<= getLastMotifColumn(); k++, pos++){
					if(pos > 0 && pos <= seqLen){
						int base = S[pos];
						_pwm[k][base]--;
						_pwm[k][0]--;
					}
				}
				it = _sites.erase(it);
				it_next++;

				/* current motif is the best => remove next motif */
			}else{
				int pos_next = it_next->pos+getFirstMotifColumn()-1;
				for(int k=getFirstMotifColumn(); k<= getLastMotifColumn(); k++, pos_next++){
					if(pos_next > 0 && pos_next <= seqLen){
						int base = S[pos_next];
						_pwm[k][base]--;
						_pwm[k][0]--;

					}
				}
				it_next = _sites.erase(it_next);
			}
		}
		if(it_next == _sites.end()) break;
	}
	updatePWM_OOL_initial(Global::pseudo);
}

void Motif::updatePval(bool Debug){
	double bestPval;
	int bestMotifNb;

	if(Debug) cerr << "update Pval: " << endl << *this;
	for(int i=0; i<2; i++){
		Kmer_Generator_Reloaded<dna_hash_t> &kmerHash = Kmer_Generator_Reloaded<dna_hash_t>::getInstance();
		kmerHash.reinitialize(getPWM(), getMotifColumns(), getPosSetSize(), Global::posSet->total[getPosSetSize()]);
		StartPosUpdater::getInstance().fill_sites_with_startPosList(kmerHash, this);
		StartPosUpdater::getInstance().fill_startPosList_with_bestSubset(this, bestPval, bestMotifNb);
		updatePWM_OOL_initial(Global::pseudo);
	}
	StartPosUpdater::getInstance().setPositionalPval(this);

	double LOG_Bonferonni = calculate_log_bonferonni(getMotifColumns(), LogTable::LOG_Neff_pwm);
	_pVal = bestPval + LOG_Bonferonni;

	if(Debug){
		cerr << "bestMotifNb: " << bestMotifNb << "\tbestPval: " << exp(_pVal) << endl;
		cerr << *this;
		//exit(-1);
	}
}

void Motif::setBindingSites(const ThresholdChecker &pValThreshold){
	Kmer_Generator_Reloaded<dna_hash_t> &kmerHash = Kmer_Generator_Reloaded<dna_hash_t>::getInstance();
	kmerHash.reinitialize(getPWM(), getMotifColumns(), getPosSetSize(), Global::posSet->total[getPosSetSize()]);
	StartPosUpdater::getInstance().fill_sites_with_possible_startPositions_oops(kmerHash, this);
	StartPosUpdater::getInstance().set_BindingSites(this, pValThreshold);
}

void Motif::updateEmpiricalPval(){
	Kmer_Generator_Reloaded<dna_hash_t> &kmerHash = Kmer_Generator_Reloaded<dna_hash_t>::getInstance();
	kmerHash.reinitialize(getPWM(), getMotifColumns(), getPosSetSize(), Global::posSet->total[getPosSetSize()]);
	StartPosUpdater::getInstance().fill_sites_with_possible_startPositions(kmerHash, 1, this);

	StartPosUpdater::getInstance().update_sites_with_empirical_Pvals(this);

	double bestPval;
	int bestMotifNb;
	StartPosUpdater::getInstance().fill_startPosList_with_bestSubset(this, bestPval, bestMotifNb);

	double LOG_Bonferonni = calculate_log_bonferonni(getMotifColumns(), LogTable::LOG_Neff_pwm);
	_pVal = bestPval + LOG_Bonferonni;

}

void Motif::buildStatePWM(double *const *const _pwm) {
	const StateLib &states = MProGlobal::getS().states;
	//	cerr << "building state pwm" << endl;
	for (int i=1; i<=_length; ++i) {
		int bestStateIndex = -1;
		double bestSum = -std::numeric_limits<double>::max();
		for (int s=0; s<(int)states.size(); ++s) {
			double sum = 0;
			for(int j=1; j<=nAlpha(Global::A); ++j) {
				//				fprintf(stderr, "states[%2d][%2d] = ", s, j);
				//				fprintf(stderr, "%f", states[s][j]);
				const double score = sqrt(exp(states[s][j]) * (_pwm[i][j]/std::max(1.0,_pwm[i][0])));
				//				fprintf(stderr, "\tadding %g for matrix cell %f/%f\n", score,_pwm[i][j],_pwm[i][0]);
				sum += score;
			}
			if (sum > bestSum) {
				bestSum = sum;
				bestStateIndex = s;
			}
			//			fprintf(stderr, "sum=%g\tbestSum=%g\tbestIndex=%d\n", sum, bestSum, bestStateIndex);
		}
		//		fprintf(stderr, "Best state %d with score %f\n", bestStateIndex, bestSum);
		for (int j=0; j<=nAlpha(Global::A); ++j) {
			_pwm[i][j] = states[bestStateIndex][j];
		}
	}
}

void Motif::buildPWM(double** pwm, double pseudo, double* bg){
	for(int i=1; i<=_length; i++){
		for(int j=1; j<= nAlpha(Global::A); j++){
			//pwm[i][j] = fast_log( ((pwm[i][j] / std::max(1.0, pwm[i][0])) + alpha*bg[j]) / (std::min(1.0, pwm[i][0]) + alpha) );
			pwm[i][j] = fast_log(static_cast<float>((pwm[i][j] + pseudo*bg[j]) / (pwm[i][0] + pseudo)));
		}
	}
}

void Motif::updatePWM_OOL_initial(double pseudo){
	double* bg = Global::negSet != NULL ? Global::negBg : Global::posBg;
	int K = getTotalSites();
	if(K==0) return;

	pseudo = getPseudocounts(pseudo, K);

	/* set whole length of motif to average of instances */
	for(int i=1; i <= PWM_LENGTH; i++){
		for(int a = 1; a <= nAlpha(Global::A); a++){_pwm[i][a] = 0;}
		_pwm[i][0] = 0; /* set to zero */
	}

	for(StartPosContainer::iterator it = _sites.begin(); it != _sites.end(); it++){
		int seq = it->seq;
		int pos = it->pos;

		updateCounts(seq, pos, _pwm);
	}

	/* set motif to optimized motif */
	buildPWM(_pwm, pseudo, bg);
}

void Motif::updatePWM_OOL(sorted_sites* sortedSites, double** pwm,
		motif_columns_type& motifColumns, int K, double pseudo){

	double* bg = Global::negSet != NULL ? Global::negBg : Global::posBg;

	pseudo = getPseudocounts(pseudo, K);

	/* set whole length of motif to average of instances */
	for(int i=1; i <= PWM_LENGTH; i++){
		for(int a = 1; a <= nAlpha(Global::A); a++){pwm[i][a] = 0;}
		pwm[i][0] = 0; /* set to zero */
	}


	int* motifsPerSequenceCount = StartPosUpdater::getInstance().getMotifsCountPerSequence();
	/* reset motifs per sequence counter */
	memset(motifsPerSequenceCount, 0, (_posSetSize+1)* sizeof(int));

	//fprintf(stderr, "update PWM with %d sequences: \n", K);
	for(int i=0, motifCounter = 0;;i++){
		int seq = sortedSites[i]->sequence;
		int pos = sortedSites[i]->startPos - motifColumns.front() + 1;

		motifsPerSequenceCount[seq]++;
		if(motifsPerSequenceCount[seq] > Global::maxMotifsPerSequence){ continue; }

		//fprintf(stderr, "%d: seq: %d pos: %d\t", i, seq, pos);
		updateCounts(seq, pos, pwm);

		if(++motifCounter == K) break;
	}

	/* set motif to optimized motif */
	buildPWM(pwm, pseudo, bg);
}

void Motif::updateCounts(int seq, int pos, double** pwm){
	int seqLen = Global::posSet->entity[seq]->n;

	int mseq = 1;

	for(int m = 0; m<mseq; m++){
		uint8_t* S = Global::posSet->entity[seq]->S[m];
		int tmpPos = pos;
		for(int j = 1; j <= _length; j++, tmpPos++){
			if(tmpPos < 1 || tmpPos > seqLen) continue;

			int base = S[tmpPos];
			//fprintf(stderr, "%c", AlphaChar(base, Global::A));
			if(base == 0) continue;
			pwm[j][base]++;
			pwm[j][0]++;
		}
		//fprintf(stderr, "\n");
	}
}

void Motif::updateCountsPWM(bool maximizeMotifLength){
	/* reset PWM */
	for(int i=1;i<=_length;i++){
		for(int j=0;j<=nAlpha(Global::A);j++)_pwm[i][j] = 0;
	}

	/* count bases */
	int start = getFirstMotifColumn();
	int	end = getLastMotifColumn();
	if(maximizeMotifLength){
		start = 1;
		end = _length;
	}

	for (StartPosContainer::const_iterator it = _sites.begin(); it !=_sites.end(); it++){
		int seq = it->seq;
		int tmpPos = it->pos;
		int n = Global::posSet->entity[seq]->n;
		if(!maximizeMotifLength) tmpPos += getFirstMotifColumn() - 1;
		int mseq = 1;

		for(int m=0; m<mseq; m++){
			uint8_t* S = Global::posSet->entity[seq]->S[m];

			int pos = tmpPos;
			for(int i=start; i<=end; i++, pos++){
				if(pos < 1 || pos > n){
					//cerr << "seq: " << it->seq << "\tstartPos: " << it->pos << endl;
					//cerr << "firstMotifColumn: " << getFirstMotifColumn() << "\tlastMotifColumn: " << getLastMotifColumn() << endl;
					//cerr << "startRegion: " << _enrichment.startRegion << "\tendRegion: "  << _enrichment.endRegion << endl;
					//cerr << pos << " not in range " << 1 << " to " << Global::posSet->entity[seq]->n << endl << endl;
					continue;
				}
				int base = S[pos];
				//if(i < 0){
				//	cerr << i << "\n";
				//	cerr << "\n";
				//}
				if(base == 0) continue;
				_pwm[i][base]++;
				_pwm[i][0]++;
			}
		}
	}
}

void Motif::InitProfMotif(char* profFile){
	for(int i=1; i<=_length; i++){
		for(int j=1; j<=nAlpha(Global::A); j++){
			_pwm[i][j] = 1.0/nAlpha(Global::A);
		}
	}

	//double pseudo = ps / (nAlpha(Global::A)-1);
	double pseudo = 0;
	FILE* fptr;
	if((fptr = fopen(profFile,"r")) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n", profFile);
		fprintf(stderr, "File does not exist!\n");
		exit(-1);
	}
	int i=0, line=1, pos = 1;
	char c='x';
	char* value = (char*)calloc(80, sizeof(char));
	while(c != EOF){
		if((c=static_cast<char>(fgetc(fptr))) == ' ' || c == '\t' || c == '\n'){}
		else break;
	}
	value[i++]=c;
	while(c != EOF){
		if((c=static_cast<char>(fgetc(fptr))) == ' ' || c == '\t' || c == '\n'){
			value[i] = '\0';
			i=0;
			if(value[0] == '-'){
				_pwm[pos][line] =  0;
			}else{
				if(line == 1) _motifColumns.push_back(pos);
				_pwm[pos][line] =  atof(value)/100.0;
			}
			if(c=='\n'){
				line++;	pos = 1;
			}else{
				pos++;
			}
			while((c=static_cast<char>(fgetc(fptr))) == ' ' || c == '\t' || c == '\n'){
				if(c=='\n'){
					line++;	pos = 1;
				}
			}
		}
		value[i++]=c;
	}


	free(value);
	int firstMotifColumn =  (_length - getMotifLength()) / 2 + 1;
	for(int i=getMotifLength(); i>0; i--){
		for(int j=1; j<=nAlpha(Global::A); j++){
			//cerr << _pwm[i][j] << "\t";
			_pwm[i+firstMotifColumn-1][j] = _pwm[i][j];
			_pwm[i+firstMotifColumn-1][0] += _pwm[i][j];
			_pwm[i][j] = 0;
		}
		//cerr << endl;
	}
	for(motif_columns_type::iterator it = getMotifColumns().begin(); it != getMotifColumns().end(); it++){
		*it += firstMotifColumn - 1;
	}


	for(i=1;i<=_length;i++){
		for(int j=1; j<= nAlpha(Global::A); j++){
			if(i>= firstMotifColumn && i <= firstMotifColumn+getMotifLength()-1){
				_pwm[i][j]=log((_pwm[i][j]+pseudo/1000)/(_pwm[i][0]+nAlpha(Global::A)*pseudo/1000));
			}else{
				_pwm[i][j]= log(1.0/nAlpha(Global::A));
			}
		}
	}

	_enrichment.max = Global::startRegion - firstMotifColumn + 1;
	_enrichment.startRegion = Global::startRegion - firstMotifColumn + 1;
	_enrichment.endRegion = Global::endRegion - firstMotifColumn + 1;

	fclose(fptr);
}

void Motif::InitStartMotif(char* startMotif, double pseudo){
	for(int i=1; i<=_length; i++){
		for(int j=1; j<=nAlpha(Global::A); j++)	_pwm[i][j] = log(1.0/nAlpha(Global::A));
	}

	const int firstMotifColumn = (_length - static_cast<int>(strlen(startMotif)) ) / 2 + 1;

	double* bg;
	bg = Global::negSet != NULL ? Global::negBg : Global::posBg;
	for(unsigned int i=firstMotifColumn;i<=firstMotifColumn + strlen(startMotif) - 1;i++){
		switch(startMotif[i-firstMotifColumn]){
		//		  		  A					C				G				T
		case 'A': _pwm[i][1]=1;		_pwm[i][2]=0;		_pwm[i][3]=0;		_pwm[i][4]=0; 		_motifColumns.push_back(i); break;
		case 'C': _pwm[i][1]=0;		_pwm[i][2]=1;		_pwm[i][3]=0;		_pwm[i][4]=0; 		_motifColumns.push_back(i); break;
		case 'G': _pwm[i][1]=0;		_pwm[i][2]=0;		_pwm[i][3]=1;		_pwm[i][4]=0; 		_motifColumns.push_back(i); break;
		case 'T': _pwm[i][1]=0;		_pwm[i][2]=0;		_pwm[i][3]=0;		_pwm[i][4]=1; 		_motifColumns.push_back(i); break;
		case 'R': _pwm[i][1]=0.5;	_pwm[i][2]=0;		_pwm[i][3]=0.5;		_pwm[i][4]=0; 		_motifColumns.push_back(i); break;
		case 'Y': _pwm[i][1]=0;		_pwm[i][2]=0.5;		_pwm[i][3]=0;		_pwm[i][4]=0.5; 	_motifColumns.push_back(i); break;
		case 'S': _pwm[i][1]=0;		_pwm[i][2]=0.5;		_pwm[i][3]=0.5;		_pwm[i][4]=0; 		_motifColumns.push_back(i); break;
		case 'W': _pwm[i][1]=0.5;	_pwm[i][2]=0;		_pwm[i][3]=0;		_pwm[i][4]=0.5;		_motifColumns.push_back(i); break;
		case 'K': _pwm[i][1]=0;		_pwm[i][2]=0;		_pwm[i][3]=0.5;		_pwm[i][4]=0.5;		_motifColumns.push_back(i); break;
		case 'M': _pwm[i][1]=0.5;	_pwm[i][2]=0.5;		_pwm[i][3]=0;		_pwm[i][4]=0; 		_motifColumns.push_back(i); break;
		case 'B': _pwm[i][1]=0;		_pwm[i][2]=1.0/3;	_pwm[i][3]=1.0/3;	_pwm[i][4]=1.0/3;	_motifColumns.push_back(i); break;
		case 'D': _pwm[i][1]=1.0/3;	_pwm[i][2]=0;		_pwm[i][3]=1.0/3;	_pwm[i][4]=1.0/3;	_motifColumns.push_back(i); break;
		case 'H': _pwm[i][1]=1.0/3;	_pwm[i][2]=1.0/3;	_pwm[i][3]=0;		_pwm[i][4]=1.0/3;	_motifColumns.push_back(i); break;
		case 'V': _pwm[i][1]=1.0/3;	_pwm[i][2]=1.0/3;	_pwm[i][3]=1.0/3;	_pwm[i][4]=0;		_motifColumns.push_back(i); break;
		case 'N': _pwm[i][1]=0.25;	_pwm[i][2]=0.25;	_pwm[i][3]=0.25;	_pwm[i][4]=0.25;	_motifColumns.push_back(i); break;
		case '.': _pwm[i][1]=0.25;	_pwm[i][2]=0.25;	_pwm[i][3]=0.25;	_pwm[i][4]=0.25;	break;
		default : cerr << "illegal character in startMotif: " << startMotif[i-1] << endl; exit(-1);
		}
	}

	/* add pseudocounts */
	double alpha = pseudo / sqrt(1000);
	for(unsigned int i=firstMotifColumn;i<=firstMotifColumn + strlen(startMotif) - 1;i++){
		for(int j=1; j<= nAlpha(Global::A); j++){
			_pwm[i][j]=log( (_pwm[i][j]+alpha*bg[j])/(1+alpha) );
		}
	}

	_enrichment.max = Global::startRegion - firstMotifColumn + 1;
	_enrichment.startRegion = Global::startRegion - firstMotifColumn + 1;
	_enrichment.endRegion = Global::endRegion - firstMotifColumn + 1;
}

void Motif::fill_startPos_with_IUPAC_matches(const std::shared_ptr<Kmer>& skr){
	/*
	 *  create an ordered list of matches
	 */
	int length = skr->getKmer()->length();
	string iupacString = skr->getKmer()->bestNucString();

	uint8_t iupac[PWM_LENGTH][3];
	uint8_t matchPos[PWM_LENGTH];
	int i=0;
	for(motif_columns_type::const_iterator it = getMotifColumns().begin(); it!= getMotifColumns().end(); it++, i++){
		const uint8_t pos = static_cast<uint8_t>(*it-getFirstMotifColumn());
		matchPos[i] = pos;
		switch(iupacString.at(pos)){
		case 'A': iupac[i][0]=1; iupac[i][1]=1;break;
		case 'C': iupac[i][0]=1; iupac[i][1]=2;break;
		case 'G': iupac[i][0]=1; iupac[i][1]=3;break;
		case 'T': iupac[i][0]=1; iupac[i][1]=4;break;
		case 'M': iupac[i][0]=2; iupac[i][1]=1;iupac[i][2]=2;break;
		case 'R': iupac[i][0]=2; iupac[i][1]=1;iupac[i][2]=3;break;
		case 'W': iupac[i][0]=2; iupac[i][1]=1;iupac[i][2]=4;break;
		case 'S': iupac[i][0]=2; iupac[i][1]=2;iupac[i][2]=3;break;
		case 'Y': iupac[i][0]=2; iupac[i][1]=2;iupac[i][2]=4;break;
		case 'K': iupac[i][0]=2; iupac[i][1]=3;iupac[i][2]=4;break;
		}
	}


	int addPositions = 1;
	int matchPositions = static_cast<int>(getMotifColumns().size());


	for(int32_t i=1; i<=(int)Global::posSet->nent; i++){

		e_type sequence = Global::posSet->entity[i];
		unsigned char* seq = sequence->S[0];

		int seqCount = 0;

		/* count hits in current sequence */
		int seqLength = sequence->n;
		//if(Global::revcomp)seqLength = sequence->n/2;
		for(int j=1; j<=seqLength; j++){
			if(Global::revcomp && j + matchPos[0] < seqLength/2 && j + matchPos[matchPositions-1] > seqLength/2) continue;
			/* check whether first match positions fit */
			for(int k=0; k<matchPositions; k++){
				if(j+matchPos[k] > seqLength)break;

				unsigned char nuc = seq[j+matchPos[k]];
				bool hit = false;
				for(int l=1; l<=iupac[k][0]; l++){
					if(nuc == iupac[k][l]){	hit = true;	break;	}
				}
				if(!hit) break;

				if(k+1==matchPositions){
					_sites.push_back(StartPos(i, static_cast<int32_t>(j-getFirstMotifColumn()+addPositions)));
					j += (length-1);   /* do not count overlapping motifs */
					seqCount++;
				}
			}
		}
	}
	//if(printDistr)fclose(fptr);
	_posSetSize = Global::posSet->nent;
}

void Motif::fill_startPos_with_seeds(const std::shared_ptr<Kmer>& skr){
	const int addPositions = 2;
	const int addSeq = 1;

	/*
	 *  create an ordered list of matches
	 */
	const MatchContainer& matches = skr->seeds;
	MatchContainer::const_iterator it_match = matches.begin();
	/* insert first element into list */
	_sites.push_back(StartPos(static_cast<int32_t>(it_match->seq + addSeq), static_cast<int32_t>(it_match->pos-getFirstMotifColumn()+addPositions)));


	StartPosContainer::iterator it = _sites.begin();
	it_match++;
	MatchContainer::const_iterator it_match_last = matches.begin();

	while(true){
		/* copy already sorted values */
		int32_t lastSeq = it_match_last->seq;
		int32_t lastPos = it_match_last->pos;
		int32_t seq = it_match->seq;
		int32_t pos = it_match->pos;
		while( (it_match != matches.end()) &&
				( (lastSeq < seq) ||
						(lastSeq == seq && lastPos < pos) ) ){
			while( it != _sites.end() &&
					( (it->seq - addSeq < seq) ||
							(it->seq - addSeq == seq && it->pos+getFirstMotifColumn()-addPositions < pos) ) ){
				it++;
			}
			_sites.insert(it, StartPos(static_cast<int32_t>(seq + addSeq), static_cast<int32_t>(pos-getFirstMotifColumn()+addPositions)));

			it_match_last++; it_match++;
			lastSeq = it_match_last->seq;
			lastPos = it_match_last->pos;
			seq = it_match->seq;
			pos = it_match->pos;
		}
		if(it_match == matches.end()) break;
		/* copy first element of new serios of sorted values */
		it = _sites.begin();
		while( it != _sites.end() &&
				( (it->seq - addSeq < seq) ||
						(it->seq - addSeq == seq && it->pos+getFirstMotifColumn()-addPositions < pos) ) ){
			it++;
		}
		_sites.insert(it, StartPos(static_cast<int32_t>(seq + addSeq), static_cast<int32_t>(pos-getFirstMotifColumn()+addPositions)));
		it_match_last++;
		it_match++;
	}

	_posSetSize = Global::posSet->nent;
}

void Motif::updateEnrichedRegion(sorted_sites* bestStartPos){
	/* update enriched region */
	resetCountStartPos();

	int* countStartPos = getCountStartPos_counts();
	int& size = getCountStartPos_size();

	int firstMotifColumn = getFirstMotifColumn();

	const int startRegion = getEnrichment().startRegion + firstMotifColumn - 1;
	const int endRegion = getEnrichment().endRegion + firstMotifColumn - 1;
	const int region_ext = static_cast<int>(round(0.5*(endRegion - startRegion + 1)));

	int i=0;

	//fprintf(stderr, "before: start: %d, end: %d, size: %d\n", getEnrichment().startRegion, getEnrichment().endRegion, size);
	//for(int i = std::max(1, startRegion - region_ext); i <= std::min(Global::posSet->max_leng, endRegion + region_ext); i++){
	//	fprintf(stderr, "%d: %d\t", i, countStartPos[i]);
	//}

	/* calculate new enriched region with all instances inside old enriched region + 50% extensions to both sides */
	while(bestStartPos[i]->score < 1e-3 || (size < 10 && bestStartPos[i]->score != 1)){
		int startPos = bestStartPos[i]->startPos;
		if(Global::revcomp && startPos > Global::posSet->max_leng / 2){
			startPos = Global::posSet->max_leng - (startPos + getMotifLength() - 1);
		}
		if(startPos >= startRegion - region_ext && startPos <= endRegion + region_ext){
			countStartPos[startPos]++;
			size++;
		}
		i++;
	}
	//fprintf(stderr, "\nsize: %d\n", size);
	//for(int i = std::max(1, startRegion - region_ext); i <= std::min(Global::posSet->max_leng, endRegion + region_ext); i++){
	//	fprintf(stderr, "%d: %d\t", i, countStartPos[i]);
	//}
	//fprintf(stderr, "\n\n");
	if(Global::usePositionalProbs){
		int maxRegion = Global::posSet->max_leng - getMotifLength() + 1;
		if(Global::revcomp){
			maxRegion = (Global::posSet->max_leng - 1) / 2 - getMotifLength() + 1;
		}
		MotifRegion motifRegion(1, maxRegion);
		region r = motifRegion.getRegion(countStartPos, size);
		r.max-= (firstMotifColumn-1) ; r.startRegion -= (firstMotifColumn-1); r.endRegion -= (firstMotifColumn-1);
		setEnrichment(r);
	}else{
		region r;
		r.set = false;
		r.max = 0;
		r.startRegion = -getFirstMotifColumn();
		r.endRegion = Global::posSet->max_leng - getMotifLength() + getFirstMotifColumn();
		setEnrichment(r);
	}

	//fprintf(stderr, "after: start: %d, end: %d\n", getEnrichment().startRegion, getEnrichment().endRegion);
}

void Motif::updateEnrichedRegion(){
	/* update enriched region */
	resetCountStartPos();

	int* countStartPos = getCountStartPos_counts();
	int& size = getCountStartPos_size();

	int firstMotifColumn = getFirstMotifColumn();
	int startRegion = getEnrichment().startRegion + firstMotifColumn - 1;
	int endRegion = getEnrichment().endRegion + firstMotifColumn - 1;
	int region_ext = static_cast<int>(round(0.5*(endRegion - startRegion + 1)));

	for(StartPosContainer::const_iterator it = getStartPosList().begin(); it != getStartPosList().end(); it++){
		int startPos = it->pos+ firstMotifColumn - 1;
		if(Global::revcomp && startPos > Global::posSet->max_leng / 2){
			startPos = Global::posSet->max_leng - (startPos + getMotifLength() - 1);
		}
		if(startPos >= startRegion - region_ext && startPos <= endRegion + region_ext){
			countStartPos[startPos]++;
			size++;
		}
	}

	if(Global::usePositionalProbs){
		int maxRegion = Global::posSet->max_leng - getMotifLength() + 1;
		if(Global::revcomp){
			maxRegion = (Global::posSet->max_leng - 1) / 2 - getMotifLength() + 1;
		}
		MotifRegion motifRegion(1, maxRegion);
		region r = motifRegion.getRegion(countStartPos, size);
		r.max-=firstMotifColumn-1; r.startRegion-=firstMotifColumn-1; r.endRegion-= firstMotifColumn-1;
		setEnrichment(r);
	}else{
		region r;
		r.set = 0;
		r.max = 0;
		r.startRegion = -getFirstMotifColumn();
		r.endRegion = Global::posSet->max_leng - getMotifLength() + getFirstMotifColumn();
		setEnrichment(r);
	}
}

Motif* Motif::getPalindrome(){
	Motif* m_palin = new Motif(*this);
	double **pwm = m_palin->getPWM();

	for(int j=getFirstMotifColumn(); j <= getLastMotifColumn(); j++){
		for(int a=1; a<=nAlpha(Global::A); a++){
			pwm[getLastMotifColumn()-j+getFirstMotifColumn()][nAlpha(Global::A)-a+1] = getPWM()[j][a];
		}
	}

	return m_palin;
}

const string Motif::getIUPACString(char* IUPAC_long) const{

	double lowThr = fast_log( 0.15f );
	double highThr = fast_log( 0.75f );

	int i, j;
	string IUPAC;
	int IUPAC_long_pos=0;
	char wildcardChar = 'N';
	for(j=getFirstMotifColumn(); j <= getLastMotifColumn(); j++){
		int count = 0;
		for(i=1; i<=nAlpha(Global::A); i++){
			//fprintf(stderr, "j: %d, i: %d, p: %f\n", j, i, exp(_pwm[j][i]));
			if(_pwm[j][i] > highThr){			// AS ueber 75% allein ausgeben
				if(IUPAC_long != NULL){
					IUPAC_long[IUPAC_long_pos++] = AlphaChar(i,Global::A);
					IUPAC_long[IUPAC_long_pos++] = ' ';
				}
				IUPAC += AlphaChar(i,Global::A);
				IUPAC += ' ';
				//fprintf(stderr, " => %s\n", IUPAC.c_str());
				count = -1;
				break;
			}else if(_pwm[j][i] > lowThr){
				count++;
			}
		}
		if(count == 0){
			if(IUPAC_long != NULL){
				IUPAC_long[IUPAC_long_pos++] = wildcardChar;
				IUPAC_long[IUPAC_long_pos++] = ' ';
			}
			IUPAC += wildcardChar;
			IUPAC += ' ';
			//fprintf(stderr, " => %s\n", IUPAC.c_str());
			continue;
		}else if(count == -1)
			continue;
		int found = 0;
		char *tmp = (char*)calloc(8, sizeof(char));
		int tmp_pos = 0;
		double max_score = -std::numeric_limits<double>::max();
		for(i=1; i<=nAlpha(Global::A); i++){
			if (_pwm[j][i]>max_score) {
				max_score = _pwm[j][i];
			}
			if(_pwm[j][i] > lowThr){			// alle AS ueber 25 % ausgeben
				if(IUPAC_long != NULL) IUPAC_long[IUPAC_long_pos++] = AlphaChar(i,Global::A);
				tmp[tmp_pos++] = AlphaChar(i, Global::A);
				if(++found < count){
					if(IUPAC_long != NULL) IUPAC_long[IUPAC_long_pos++] = '/';
				}
			}
		}
		tmp[tmp_pos] = '\0';
		IUPAC += getIUPAC_char(tmp, Global::A);
		IUPAC += ' ';
		//fprintf(stderr, " => %s\n", IUPAC.c_str());
		free(tmp);
		if(IUPAC_long != NULL) IUPAC_long[IUPAC_long_pos++] = ' ';
	}
	if(IUPAC_long != NULL) IUPAC_long[IUPAC_long_pos] = '\0';

	motif_columns_type mcols = getMotifColumns();
	int col = getFirstMotifColumn();
	for (size_t i=0; i<IUPAC.length(); ++i) {
		if (IUPAC[i] == ' ') {
			continue;
		}
		bool found = false;
		for (motif_columns_type::const_iterator c=mcols.begin(); !found && c!=mcols.end(); ++c) {
			if (*c==(int)col) {
				found = true;
			}
		}
		if (!found) {
			IUPAC[i] = '.';
		}
		++col;
	}

	return IUPAC;
}

void Motif::makeReverseComplement() {
	//  assert(!Global::aa); /* doesn't make sense for proteins */
	assert(true);
	for (int i=getFirstMotifColumn(), j=getLastMotifColumn(); i<=j; ++i,--j) {
		std::swap(_pwm[i][AlphaCode('A', Global::A)], _pwm[j][AlphaCode('T', Global::A)]);
		std::swap(_pwm[i][AlphaCode('C', Global::A)], _pwm[j][AlphaCode('G', Global::A)]);
		if (i<j) {
			std::swap(_pwm[i][AlphaCode('G', Global::A)], _pwm[j][AlphaCode('C', Global::A)]);
			std::swap(_pwm[i][AlphaCode('T', Global::A)], _pwm[j][AlphaCode('A', Global::A)]);
		}
	}
}

void Motif::printPWM(ostream &os, double** pwm) const{
	os << "\t\t";
	for(int i=0;i<getLastMotifColumn()-getFirstMotifColumn()+1;i++)
		os << std::setprecision(2) << std::setw(6) << i+1 << "\t";
	os << endl << "\t-------";
	for(int i=getFirstMotifColumn();i<=getLastMotifColumn();i++)os << "--------";
	for(int i=1;i<=nAlpha(Global::A);i++){
		os << endl << "\t" << AlphaChar(i,Global::A) << "\t";
		os << std::resetiosflags( std::ios::scientific );
		os << std::fixed;
		motif_columns_type::const_iterator it = getMotifColumns().begin();
		while(true){
			int column = *it;
			os << std::setw(6) << (int)(exp(pwm[*it][i]) * 10000)/100.00 << "\t";
			it++;
			if(it == getMotifColumns().end())break;
			for(int col = column+1; col < *it; col++) os << std::setw(6) << "-" << "\t";
		}
	}
	os << endl << endl;
}

void Motif::printLogPWM(ostream &os, double** pwm, motif_columns_type& motifColumns, bool log) const{
	os << "\t\t";
	for(int i=0;i<motifColumns.back()-motifColumns.front()+1;i++)
		os << std::setprecision(2) << std::setw(6) << i+1 << "\t";
	os << endl << "\t-------";
	for(int i=motifColumns.front();i<=motifColumns.back();i++)os << "--------";
	for(int i=1;i<=nAlpha(Global::A);i++){
		os << endl << "\t" << AlphaChar(i,Global::A) << "\t";
		os << std::resetiosflags( std::ios::scientific );
		os << std::fixed;
		motif_columns_type::const_iterator it = motifColumns.begin();
		while(true){
			int column = *it;

			if(log)  os << std::setw(6) << (int)(exp(pwm[*it][i])*10000)/100.00 << "\t";
			else os << std::setw(6) << (int)(pwm[*it][i]* 10000)/100.00 << "\t";

			it++;
			if(it == motifColumns.end())break;
			for(int col = column+1; col < *it; col++) os << std::setw(6) << "-" << "\t";
		}
	}
	os << endl << endl;
}

void Motif::printFullPWM(ostream &os) const{
	os << endl << "\t\t";
	for(int i=0;i<getLength();i++)
		os << std::setprecision(2) << std::setw(6) << i+1 << "\t";
	os << endl << "\t-------";
	for(int i=0;i<=_length;i++)os << "--------";
	for(int i=1;i<=nAlpha(Global::A);i++){
		os << endl << "\t" << AlphaChar(i,Global::A) << "\t";
		os << std::resetiosflags( std::ios::scientific );
		os << std::fixed;
		for(int column=1; column<=getLength(); column++){
			os << std::setw(6) << (int)(exp(getPWM()[column][i]) * 10000)/100.00 << "\t";
		}
	}
	os << endl << endl;
}

ostream& operator<< (ostream &os, const Motif& M){
	std::iostream::fmtflags flags = os.flags();
	std::streamsize prec = os.precision();

	if(M.getTotalSites() > 0){
		os << endl << "\t======="; for(int i=M.getFirstMotifColumn();i<=M.getLastMotifColumn();i++)os << "========";
		os << endl;

		double pVal_log10 = M.getPval() / LogTable::LOG_10;

		os << "\t" << M.getIUPACString() << ":" << endl <<
				"\t  " << M.getTotalSites() << " sites in " << M.getPosSetSize() << " sequences" << endl;
		//if(M.getTotalBindingSites() > 0) os << "\t Binding Sites: " << M.getTotalBindingSites();
		if(!Global::multipleOccurrence){
			os << std::setprecision(2) <<
					"\t => occurrence: " <<  M.getTotalSites() *100.0/ M.getPosSetSize() << "%" << std::scientific << endl;
		}
		//os << "\t  log(E-Value): " << M._pVal << "\t";

		if (pow(10, pVal_log10) == 0) os << "\t  E-Value: 1e" << (int)pVal_log10 << endl;
		else os << "\t  E-Value: " << pow(10, pVal_log10) << "\t";

		os << "\t pos Pval: " << pow(10, M.getPosPval()) << endl;
		int offset = Global::posSet->max_leng - Global::downstream;
		if(Global::usePositionalProbs ){
			/* test whether a significant region could be found */
			if(M._enrichment.set == 0){
				os << "\t  max: - \tRegion: -/-"<< endl;
			}else{
				os << "\t  max: " << M._enrichment.max + M.getFirstMotifColumn() -1 -offset <<
						"\tRegion: " << M._enrichment.startRegion + M.getFirstMotifColumn() -1 -offset << "/"<<
						M._enrichment.endRegion + M.getFirstMotifColumn() -1 -offset << endl;
			}
		}

		//    os << endl;
		//    for(StartPosContainer::const_iterator it = M._sites.begin(); it != M._sites.end(); it++){
		//      //int offset = Global::posSet->entity[it->seq]->n - Global::downstream;
		//      os << it->seq << "/" << it->pos /*+ M.getFirstMotifColumn() -1 - offset*/ << " ";
		//    }
		//    os << endl;


		//    cerr << endl;
		//    for(StartPosContainer::const_iterator it = M._sites.begin(); it != M._sites.end(); it++){
		//      //if(it->seq != 1)continue;
		//      for(int m=0; m<1; m++){
		//        //for(int m=0; m<Global::posSet->entity[it->seq]->mseq; m++){
		//        fprintf(stderr, "%d/%d (%d)\t", it->seq, it->pos, Global::posSet->entity[it->seq]->n);
		//        for(int j=0; j<M.getMotifLength(); j++){
		//          fprintf(stderr, "%c", AlphaChar(Global::posSet->entity[it->seq]->S[m][it->pos+M.getFirstMotifColumn()+j-1], Global::A));
		//        }
		//        cerr << endl;
		//      }
		//      //cerr << endl;
		//    }
		//    cerr << endl;
	}

	os << endl << "\tMODEL:" << endl;
	M.printPWM(os, M.getPWM());

	os.flags(flags);
	os << std::setprecision(static_cast<int>(prec));

	return os;
}

void Motif::calculateInterpolatedMarkovModelProbs(	bool lastCondsPseudoCounts, int minorder /* = 0 */ ){

	int order = minorder;
	if( order < 1 ){
		if( order <= _order ){

			motif_columns_type::iterator iter;
			// calculate 0th-order counts numbers
			// total counts differ between columns in case sequences contain Ns
			double* ncounts = ( double* )calloc( _length+1, sizeof( double ) );
			for( iter=_motifColumns.begin(); iter != _motifColumns.end(); iter++
			){
				for( int i=1; i <= nAlpha( Global::A ); i++ ){
					ncounts[*iter] += _counts[*iter][i];
				}
			}

			// calculate 0th-order probabilities
			for( iter=_motifColumns.begin(); iter != _motifColumns.end(); iter++
			){
				for( int k=1; k <= nAlpha( Global::A ); k++ ) {
					_pwm[*iter][k] = ( _counts[*iter][k] + _alpha[0] *
							Global::posBg[k] ) / ( ncounts[*iter] +
									_alpha[0] );
					_conds[*iter][k] = _pwm[*iter][k];
				}
			}
			free( ncounts );

			order++;
		}
	}

	// calculate higher-order probabilities
	uint8_t* kmer = ( uint8_t* )malloc( ( _order+1 )*sizeof( uint8_t ) );
	for( ; order <= _order; order++ ){
		calculateHoInterpolatedMarkovModelProbs( kmer, 0, order,
				lastCondsPseudoCounts );
	}
	free( kmer );
}

void Motif::calculateMarkovModelProbs(){

	int order = 0;
	if( order <= _order ){

		motif_columns_type::iterator iter;

		// calculate 0th-order counts
		// total counts differ between columns in case sequences contain Ns
		double* ncounts = ( double* )calloc( _length+1, sizeof( double ) );
		for( iter=_motifColumns.begin(); iter != _motifColumns.end(); iter++ ){
			for( int i=1; i <= nAlpha( Global::A ); i++ ){
				ncounts[*iter] += _counts[*iter][i];
			}
		}
		// calculate 0th-order probabilities
		for( iter=_motifColumns.begin(); iter != _motifColumns.end(); iter++ ){
			for( int k=1; k <= nAlpha( Global::A ); k++ ) {
				_pwm[*iter][k] =
						( _counts[*iter][k] + _alpha[0]*Global::posBg[k] ) /
						( ncounts[*iter] + _alpha[0] );
				_conds[*iter][k] = _pwm[*iter][k];
			}
		}
		free( ncounts );

		// calculate higher-order probabilities
		uint8_t* kmer = ( uint8_t* )malloc( ( _order+1 )*sizeof( uint8_t ) );
		for( order++; order <= _order; order++ ){
			calculateHoMarkovModelProbs( kmer, 0, order );
		}
		free( kmer );
	}
}

void Motif::initHoMotifWithBindingSites( char* bindingSiteFile ){

	int k, l, L;
	char character;
	char* line;
	unsigned char* sequence;

	FILE* fp;
	list<int>::iterator iter;

	int b = -1;
	int bindingSiteLength = _length - Global::addColumns.at( 0 ) -
			Global::addColumns.at( 1 );

	float weight;
	std::vector<float> weights;

	if( Global::bindingSiteIntsFile != NULL ){
		std::ifstream fs( Global::bindingSiteIntsFile, std::ios_base::in );
		if( fs.fail() ){
			fprintf( stderr, "Cannot open bindingSiteIntsFile %s\n",
					Global::bindingSiteIntsFile );
			exit( -1 );
		}

		while( fs >> weight ){
			if( weight < 0.0f ){
				fprintf( stderr, "Use non-negative binding site weights\n" );
				exit( -1 );
			}
			weights.push_back( weight );
		}
		fs.close();

		// calculate weights
		if( Global::rankWeighting ){

			// calculate rank-based weights
			calculateWeights( weights, Global::bindingSiteBackgroundRank, true
			);
		} else{

			// calculate intensity-based weights
			calculateWeights( weights, Global::bindingSiteBackgroundIntensity,
					false );
		}
	}

	if( ( fp = fopen( bindingSiteFile, "r" ) ) == NULL ){
		fprintf( stderr, "Cannot open bindingSiteFile %s\n", bindingSiteFile
		);
		exit(-1);
	}

	line = ( char* )malloc( bindingSiteLength+2 * sizeof( char ) ); // +2: \n and \0
	sequence = ( unsigned char* )malloc( bindingSiteLength+1 * sizeof( uint8_t ) );

	int sequenceCounter = 0;

	/* calculate counts from sequences */

	/* read in binding sites line-by-line
	 * and skip blank lines */
	while( fgets( line, bindingSiteLength+2, fp ) != NULL ){ // +2: \n and \0

		if( line[0] == '\n' ){
			continue; // skip blank lines
		} else{
			sequenceCounter++; // sequence counter

			if( Global::bindingSiteIntsFile != NULL ){
				if( sequenceCounter <= static_cast<int>( weights.size() ) ){
					weight = weights.at( sequenceCounter-1 );
				} else{
					fprintf( stderr, "The number of binding sites and binding "
							"site weights differs\n" );
					exit( -1 );
				}
			} else{
				weight = 1.0f;
			}
		}

		/* convert to upper-case letter
		 * and copy to sequence without \n and \0 */
		for( L=0; ( character = line[L] ) != '\n'; ++L ){
			sequence[L] = AlphaCode( ( islower(character) )? toupper(character)
					: character, Global::A );
		}

		if( L != bindingSiteLength ){
			fprintf( stderr, "Binding site sequence lengths must not "
					"differ\n" );
			exit(-1);
		}

		if( b < 0 ){
			// set motif indices in data structure (starting at 1)
			for( int l=1; l <= _length; l++ ){
				_motifColumns.push_back( l ); // set motif indices
			}
			// set binding site offsets
			b = Global::addColumns.at( 0 ) + 1;
		}

		for( l=0; l < L; l++ ){ // across positions
			for( k=0; k <= _order && l+k < L && sequence[l+k]; k++ ){
				/* sequence[l+k]
				 * checks whether nucleotide equals N
				 * represented by 0 */
				_counts[b+l+k][sub2ind( sequence, l, k )] += weight;
				if( k > 0 ){
					// k-mer not followed by N character
					_countsx[b+l+k-1][sub2ind( sequence, l, k-1 )] += weight;
				}
			}
		}
	}
	if( sequenceCounter < static_cast<int>( weights.size() ) ){
		fprintf( stderr, "The number of binding sites and binding site weights "
				"differs\n" );
		exit( -1 );
	}

	if( Global::msq ){
		setSpecificityFactor( std::min( static_cast<float>( sequenceCounter ) /
				static_cast<float>( Global::posSet->nent ),
				Global::qmax ) );
	}

	if( Global::interpolate ){
		// calculate interpolated Markov model probabilities from counts
		calculateInterpolatedMarkovModelProbs( false, false );
	} else{
		// calculate Markov model probabilities from counts
		calculateMarkovModelProbs();
	}


	if( Global::verbose ){
		printInterpolatedMarkovModel( *this, true );
	}

	free( line );
	free( sequence );
}

void Motif::initHoMotifWithInterpolatedMarkovModel( char* baseFileName ){

	// set motif indices in data structure (starting at 1)
	for( int l=1; l <= _length; l++ ){
		_motifColumns.push_back( l ); // set motif indices
	}

	// set binding site offsets
	int b = Global::addColumns.at( 0 ) + 1;

	FILE* fp;
	std::stringstream str;

	// read in conditional probababilities

	str << baseFileName << ".conds";

	if( ( fp = fopen( str.str().c_str(), "r" ) ) == NULL ){
		fprintf( stderr, "Cannot open markovModelFile %s\n", str.str().c_str()
		);
		exit(-1);
	}

	float value; // probability

	int pos = b;
	int field = 0;

	while( fscanf( fp, "%e", &value ) != EOF ){
		++field;
		if( field == _fields ){
			++pos;
			field = 1;
		}
		_conds[pos][field] = value;
	}

	fclose( fp );

	// read in probababilities

	str.str( "" );
	str << baseFileName << ".probs";

	if( ( fp = fopen( str.str().c_str(), "r" ) ) == NULL ){
		fprintf( stderr, "Can't open file %s\n", str.str().c_str() );
		exit(-1);
	}

	pos = b;
	field = 0;

	while( fscanf( fp, "%e", &value ) != EOF ){
		++field;
		if( field == _fields ){
			++pos;
			field = 1;
		}
		_pwm[pos][field] = value;
	}

	fclose( fp );


	if( Global::verbose ){
		printInterpolatedMarkovModel( *this, true );
	}
}

void Motif::initHoMotifWithStartPos(){

	if( getTotalSites() < 1 ){
		fprintf( stderr, "No binding site instances to initialize "
				"higher-order model.\n" );
		exit( -1 );
	}

	_order = std::min( _motifColumns.back()-_motifColumns.front(), _order );

	if( _order < _motifColumns.back()-_motifColumns.front() ){

		_alpha.resize( _order+1 );
	}


	_coeffs = ( int* )calloc( _order+2, sizeof( int ) );
	_offsets = ( int* )calloc( _order+2, sizeof( int ) );

	int k;
	for( k=0; k < _order+2; k++ ){
		_coeffs[k] = static_cast<int>( pow( nAlpha( Global::A ), k ) );
		_offsets[k] = k ? _offsets[k-1]+_coeffs[k] : _coeffs[k];
	}
	_fields = _offsets[k-1];

	for( int i=0; i <= _length; i++ ){
		free( _pwm[i] );
	}
	free( _pwm );

	int L = getMotifLength();
	int offset = Global::addColumns.at( 0 ) + 1;

	_length = Global::addColumns.at( 0 ) + L + Global::addColumns.at( 1 );

	_pwm = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
	for( int i=0; i <= _length; i++ ){
		_pwm[i] = ( double* )calloc( _fields, sizeof( double ) );
	}

	_counts = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
	for( int i=0; i <= _length; i++ ){
		_counts[i] = ( double* )calloc( _fields, sizeof( double ) );
	}
	_countsx = ( double** )malloc( ( _length )*sizeof( double* ) );
	for( int i=0; i < _length; i++ ){
		_countsx[i] = ( double* )calloc( _offsets[_order], sizeof( double ) );
	}

	_conds = ( double** )malloc( ( _length+1 )*sizeof( double* ) );
	for( int i=0; i <= _length; i++ ){
		_conds[i] = ( double* )calloc( _fields, sizeof( double ) );
	}
	_lastConds = NULL;

	int b = _motifColumns.front();

	_motifColumns.clear();
	for( int i=1; i <= _length; i++ ){
		_motifColumns.push_back( i );
	}

	float weight = 1.0f;
	StartPosContainer::iterator iter;
	for( iter=_sites.begin(); iter != _sites.end(); iter++ ){

		uint8_t* sequence = Global::posSet->entity[iter->seq]->S[0];

		if( Global::initInts ){
			weight = Global::posSet->entity[iter->seq]->weight;
		}

		int startPos = iter->pos + b - 1;
		int endPos = startPos + L;

		for( int l=0; l < L; ++l, ++startPos ){
			for( k=0; k <= _order && startPos+k < endPos &&
			sequence[startPos+k]; k++ ){
				/* sequence[startPos+k]
				 * checks whether nucleotide equals N
				 * represented by 0 */
				_counts[offset+l+k][sub2ind( sequence, startPos, k )] += weight;
				if( k > 0 ){
					// k-mer not followed by N character
					_countsx[offset+l+k-1][sub2ind( sequence, startPos, k-1 )]
					                       += weight;
				}
			}
		}
	}

	if( Global::interpolate ){
		// calculate interpolated Markov model probabilities from counts
		calculateInterpolatedMarkovModelProbs( false );
	} else{
		// calculate Markov model probabilities from counts
		calculateMarkovModelProbs();
	}


	if( Global::verbose ){
		printInterpolatedMarkovModel( *this, true );
	}
}

void Motif::setPseudoCounts( std::vector<float> a ){
	_alpha = a;
}

/*
 * sets the percentage of positive sequences containing a binding site
 * instance as the model-specific specificity factor used in EM calculations
 * and recalculates pseudocounts factor
 */
void Motif::setSpecificityFactor( float q ){

	_q = q;
}

// calculateKmerProbs( kmer[order+1], 0, order )
double Motif::calculateHoLogPrior( unsigned char* kmer, int pos, int order ){

	double prior = 0.0;

	for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
		kmer[pos] = k; // next kmer base
		if( pos == ( order-1 ) ){
			// calculate prior
			prior += calculateHoLogPrior( kmer, order );
		} else{ // add another base
			prior += calculateHoLogPrior( kmer, pos+1, order );
		}
	}

	return prior;
}

double Motif::calculateHoLogPrior( unsigned char* kmer, int order ){

	int i, ii;
	int modelPos;
	double prior = 0.0;

	list<int>::iterator iter;
	for( iter=_motifColumns.begin(); iter != _motifColumns.end(); iter++ ){

		modelPos = *iter - *_motifColumns.begin(); // 0, 1, ...

		if( modelPos < order ){ // margins
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// shorter kmer indices
				i = sub2ind( kmer+(order-modelPos), modelPos );
				// longer kmer indices
				ii = sub2ind( kmer, order );
				// calculate prior
				prior += _conds[*iter][i] * log( _conds[*iter][ii]);
			}
		}
		else{ // remaining model positions
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// shorter kmer indices
				i = sub2ind( kmer+1, order-1 );
				// longer kmer indices
				ii = sub2ind( kmer, order );
				// calculate prior
				prior += _conds[*iter][i] * log( _conds[*iter][ii] );
			}
		}
	}

	return prior;
}

void Motif::calculateHoInterpolatedMarkovModelProbs( unsigned char* kmer,
		int pos, int order, bool lastCondsPseudoCounts ){

	for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
		kmer[pos] = k; // next kmer base
		if( pos == ( order-1 ) ){
			// calculate probabilities
			calculateHoInterpolatedMarkovModelProbs( kmer, order,
					lastCondsPseudoCounts );
		} else{ // add another base
			calculateHoInterpolatedMarkovModelProbs( kmer, pos+1, order,
					lastCondsPseudoCounts );
		}
	}
}

void Motif::calculateHoInterpolatedMarkovModelProbs( unsigned char* kmer,
		int order, bool lastCondsPseudoCounts ){

	int i, ii, p;
	int modelPos;

	list<int>::iterator iter;
	for( iter=_motifColumns.begin(); iter != _motifColumns.end(); iter++ ){

		modelPos = *iter - *_motifColumns.begin(); // 0, 1, ...

		if( modelPos < order ){ // margins
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// shorter kmer indices
				i = sub2ind( kmer+(order-modelPos), modelPos );
				// longer kmer indices
				ii = sub2ind( kmer, order );
				// defaults to lower-order conditional probabilities
				_conds[*iter][ii] = _conds[*iter][i];
				// defaults to lower-order probabilities
				_pwm[*iter][ii] = _pwm[*iter][i];
			}
		}
		else{ // remaining model positions
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// shorter kmer indices
				i = sub2ind( kmer, order-1 );
				// longer kmer indices
				ii = sub2ind( kmer, order );
				// pseudo-counts kmer indices
				p = sub2ind( kmer+1, order-1 );
				// calculate higher-order conditional probabilities
				if( lastCondsPseudoCounts ){
					// using last iteration's pseudo-counts
					_conds[*iter][ii] =
							( _counts[*iter][ii] + _alpha[order]*_lastConds[*iter][p] )
							/ ( _countsx[*iter-1][i] + _alpha[order] );
				} else{
					_conds[*iter][ii] =
							( _counts[*iter][ii] + _alpha[order]*_conds[*iter][p] )
							/ ( _countsx[*iter-1][i] + _alpha[order] );
				}
				// calculate higher-order probabilities
				_pwm[*iter][ii] = _conds[*iter][ii] * _pwm[*iter-1][i];
			}
		}
	}
}

// calculateKmerProbs( kmer[order+1], 0, order, 1.0 )
void Motif::calculateHoMarkovModelProbs( unsigned char* kmer, int pos, int order
){

	for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
		kmer[pos] = k; // next kmer base
		if( pos == ( order-1 ) ){
			// calculate probabilities
			calculateHoMarkovModelProbs( kmer, order );
		} else{ // add another base
			calculateHoMarkovModelProbs( kmer, pos+1, order );
		}
	}
}

void Motif::calculateHoMarkovModelProbs( unsigned char* kmer, int order ){

	int i, ii;
	int modelPos;

	list<int>::iterator iter;
	for( iter=_motifColumns.begin(); iter != _motifColumns.end(); iter++ ){

		modelPos = *iter - *_motifColumns.begin(); // 0, 1, ...

		if( modelPos < order ){ // margins
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// shorter kmer indices
				i = sub2ind( kmer+(order-modelPos), modelPos );
				// longer kmer indices
				ii = sub2ind( kmer, order );
				// lower-order conditional probabilities
				_conds[*iter][ii] = _conds[*iter][i];
				// lower-order probabilities
				_pwm[*iter][ii] = _pwm[*iter][i];
			}
		}
		else{ // remaining model positions
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// shorter kmer indices
				i = sub2ind( kmer, order-1 );
				// longer kmer indices
				ii = sub2ind( kmer, order );
				// calculate higher-order conditional probabilities
				_conds[*iter][ii] =
						( _counts[*iter][ii] + _alpha[order]*Global::posBg[k] )
						/ ( _countsx[*iter-1][i] + _alpha[order] );
				// calculate higher-order probabilities
				_pwm[*iter][ii] = _conds[*iter][ii] * _pwm[*iter-1][i];
			}
		}
	}
}
void Motif::resetMotifParams(int order, double alpha){
	_Qfunc_p1 = 0.0;
	_Qfunc_p2 = 0.0;
	_Qfunc_p3 = 0.0;
	_Qfunc_p4 = 0.0;

	_Qfunc_grad = 0.0;
	_Grad_p1    = _length * pow(4,order) * double(digamma(alpha + 4));
	_Grad_p2    = 0.0;
	_Grad_p3    = 0.0;

	_zeroOrder   = 0.0;
	_pseudoOrder = 0.0;

}

void Motif::updateConds(unsigned char* kmer, int order, double alpha,int pos,double* counts_all){
	if(order == 0){
		int y;
		list<int>::iterator j;

		for(j = this->getMotifColumns().begin(); j != this->getMotifColumns().end(); j++){
			for( unsigned char k_0=1; k_0 <= nAlpha( Global::A ); k_0++ ){
				kmer[0] = k_0;
				// longer kmer indices
				y = sub2ind( kmer,order );
				_conds[*j][y] = (_counts[*j][y] + alpha * Global::posBg[k_0])/(counts_all[*j]  + alpha);

			}//end for loop k
		}// end for loop j=1..W-1
	}
	else{
		for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
			kmer[pos] = k; // next kmer base
			if( pos == ( order-1 )){
				// calculate score
				updateConds(kmer,order, alpha);
			} else{
				// add another base
				updateConds(kmer,order, alpha,pos+1,counts_all);
			}
		}
	}
}
void Motif::updateConds( unsigned char* kmer, int order, double alpha){

	int y,y_prime,p, modelPos;
	list<int>::iterator j;
	double** conds = _conds;

	for(j = this->getMotifColumns().begin(); j != this->getMotifColumns().end(); j++){

		//modelPos = *j - this->getMotifColumns().begin(); //0,1,...
		modelPos = *j - 1; //0,1,...
		for( unsigned char k_0=1; k_0 <= nAlpha( Global::A ); k_0++ ){
			kmer[order] = k_0;
			if( modelPos < order ){ // margins
				// shorter kmer indices
				y_prime = sub2ind( kmer+(order-modelPos), modelPos);
				// longer kmer indices
				y = sub2ind( kmer,order );
				// defaults to lower-order contidional probabilities
				_conds[*j][y] = conds[*j][y_prime];
			}// end if
			else{ // remaining model positions
				// shorter kmer indices
				y_prime = sub2ind( kmer,order-1 );
				// longer kmer indices
				y = sub2ind( kmer,order );
				//pseudo-counts kmer indices
				p = sub2ind(kmer+1,order-1);
				_conds[*j][y] = (_counts[*j][y] + alpha * conds[*j][p])/(_counts[*j-1][y_prime] + alpha);
			}//end else
		}//end for loop k
	}// end for loop j=1..W-1
}


