#ifndef _BB_INL_H
#define _BB_INL_H

#include "branch_and_bound.h"
#include "NullModel.h"
#include "UngappedKmer.h"
#include "memoryPool/pool_alloc.h"
#include <iostream>


template<class HASH>
Kmer_Generator_Reloaded<HASH>::Kmer_Generator_Reloaded() : bgLog((Global::negSet==NULL) ? Global::posBg_log : Global::negBg_log) {
	ID_SHIFT = 2;
	CHR_MASK = 3;
	PWM = new double*[PWM_LENGTH+1];
	PWM_lin = new double*[PWM_LENGTH+1];

	debug = false;

	for (int i = 0; i <= PWM_LENGTH; ++i) {
		PWM[i] = new double[nAlpha(Global::A)+1];
		PWM_lin[i] = new double[nAlpha(Global::A)+1];
	}

}



template<class HASH>
void Kmer_Generator_Reloaded<HASH>::reinitialize_common() {
	gaps.clear();
	pValueHash[0].clear();
	pValueHash[1].clear();
	probabilityHash.clear();
	scoreIdVector[0].clear();
	scoreIdVector[1].clear();
	scoreMapArray.clear();
}



template<class HASH>
inline int Kmer_Generator_Reloaded<HASH>::get_Kmer_score(unsigned char* kmer){
	double score = 0;
	for(int i=1; i<=totalLength; i++){
		score += PWM[i][kmer[i-1]];
	}
	/* - 1 for rounding errors */
	return (int)(score*1e6 - 1);
}



template<class HASH>
inline double Kmer_Generator_Reloaded<HASH>::calculate_kmer_probability(const id_type kmer_id){
	/* Generate 'expanded' kmer, i.e., a string with explicit gaps,
	 * then use NullModel for calculating its probability.
	 */
	static unsigned char expanded_kmer[100];
	id_type id = kmer_id;
	int totalPos = 0;
	for (int matchPos=0; matchPos<length; ++matchPos) {
		/* indices in id are shifted by -1 w.r.t. the alphabet */
		expanded_kmer[totalPos++] = static_cast<uint8_t>((id & CHR_MASK) + 1);
		id >>= ID_SHIFT;
		for (int i=0; i<gaps[matchPos+gapOffset]; ++i) {
			expanded_kmer[totalPos++] = 0;
		}
	}

	/* calculate probability by using calculated PWM dependent
	 * conditionals for crossing the split
	 */
	double prob;
	if (part==0) {
		prob = NullModel::getProbability(expanded_kmer, totalPos);
		if(fillProbabilityHash) {
			prefixList.push_back(pwm_kmer_t(UngappedKmer(expanded_kmer, totalPos)));
		}
		//fprintf(stderr, "\t%e\n", prob);
	} else {
		prob = 1;
		for (int stop = 0; stop < totalPos; ++stop) {
			const int start = std::max(0, stop - NullModel::getOrder());
			const int len = stop-start+1;
			const UngappedKmer k(expanded_kmer+start, len);
			//cerr << k;
			if (len <= NullModel::getOrder()) {
				double rightProb = getRightKmerConditional(k);
				prob *= rightProb;
			} else {
				prob *= NullModel::getConditional(k);
			}
		}
		//fprintf(stderr, "\t%e\n", prob);
		//exit(-1);
	}

	return prob;
}


template<class HASH>
inline double Kmer_Generator_Reloaded<HASH>::get_kmer_probability_matched_only(const std::string &kmer) {
	static unsigned char k[100];
	assert((int)kmer.length() == totalLength);
	for (size_t i=0; i<kmer.length(); ++i) {
		k[i] = AlphaCode(kmer[i], Global::A);
	}
	return get_kmer_probability_matched_only(k);
}





template<class HASH>
void Kmer_Generator_Reloaded<HASH>::reinitialize(const AbstractKmer &kmer,
		const int splitAfter, const double threshold_left, const double threshold_right) {
	reinitialize_common();

	MAX_KMER_SIZE = splitAfter;
	hashNumber = 2;
	totalLength = kmer.numMatches();
	if(totalLength <= MAX_KMER_SIZE) hashNumber = 1;

	/* create gaps and PWM data structures */
	gaps.resize(totalLength);
	for (int kpos=0; kpos<totalLength-1; ++kpos) {
		gaps[kpos]= kmer.gapsAfter(kpos);
	}
	gaps[totalLength-1] = 0;
	int i=1;
	for (int kpos=0; kpos<totalLength; ++kpos) {
		int kchar = kmer.charAt(kpos);
		double lin_sum = 0;
		for (int j = 1; j <= nAlpha(Global::A); ++j) {
			PWM[i][j] = MProGlobal::getS().states[kchar][j] - bgLog[j];
			PWM_lin[i][j] = exp(MProGlobal::getS().states[kchar][j]);
			lin_sum += PWM_lin[i][j];
		}
		assert(std::abs(lin_sum - 1) < 1e-3);
		++i;
	}

	//calculateRightProbabilities();

	/* do not split kmer in several parts */
	if(totalLength <= MAX_KMER_SIZE){
		//pValueHash[0].rehash(pow(4,totalLength+2));
		part = 0;
		fillProbabilityHash = false;
		T = (int)(threshold_left*1e6 - 1);
		gapOffset = 0;
		set_kmer_to_Pvalue_hash(1, totalLength);

	/* split kmer in several parts */
	}else{
		T = (int)(threshold_left*1e6 - 1);
		prefixList.clear();
		part = 0;
		fillProbabilityHash = true;
		gapOffset = 0;

		set_kmer_to_Pvalue_hash(1, MAX_KMER_SIZE);
		scoreMapArray.clear();

		T = (int)(threshold_right*1e6 - 1);
		fillProbabilityHash = false;
		part = 1;

		recalculatePrefixProbabilities(prefixList);

		int start = MAX_KMER_SIZE+1;
		int stop = totalLength;
		gapOffset = MAX_KMER_SIZE;
		set_kmer_to_Pvalue_hash(start, stop);
		scoreMapArray.clear();
	}
}


template<class HASH>
void Kmer_Generator_Reloaded<HASH>::reinitialize(double const* const* pwm, const motif_columns_type& motif_columns, int pos_set_size, int total_positions) {
	reinitialize_common();
	posSetSize = pos_set_size;
	totalPositions = total_positions;

	totalLength = static_cast<int>(motif_columns.size());
	MAX_KMER_SIZE = 8;
	if(totalLength > 17){
		MAX_KMER_SIZE = totalLength / 2;
	}

	/* set score threshold which has to be satisfied by the kmers expected score if kmer matches a random matrix */
	hashNumber = 2;
	if(totalLength <= MAX_KMER_SIZE) hashNumber = 1;

	/* fill gaps and PWM data structures */
	gaps.resize(totalLength);
	std::vector<int> _motifColumns(totalLength);
	motif_columns_type::const_iterator it = motif_columns.begin();
	_motifColumns[0] = *it;
	it++;
	for(int i=1; it != motif_columns.end(); it++, i++){
		_motifColumns[i] = *it;
		gaps[i-1] = _motifColumns[i] - _motifColumns[i-1] - 1;
	}
	gaps[totalLength-1] = 0; /* 'zero' gap after last residue */
	for (int i=1; i<=totalLength; ++i) {
		double lin_sum = 0;
		for (int j = 1; j <= nAlpha(Global::A); ++j) {
			PWM[i][j] = pwm[_motifColumns[i-1]][j] - bgLog[j];
			PWM_lin[i][j] = exp(pwm[_motifColumns[i-1]][j]);
			lin_sum += PWM_lin[i][j];
		}
		if(std::abs(lin_sum-1) >= 1e-3){
			for(int j= 1; j<= nAlpha(Global::A); ++j){
				fprintf(stderr, "%f\t", PWM_lin[i][j]);
			}
			fprintf(stderr, "\nlin_sum: %f\n", lin_sum);
		}
		assert(std::abs(lin_sum - 1) < 1e-3);
	}

	//calculateRightProbabilities();

	/* do not split kmer in several parts */
	if(totalLength <= MAX_KMER_SIZE){
		part = 0;
		fillProbabilityHash = false;
		gapOffset = 0;
		T = 0;
		set_kmer_to_Pvalue_hash(1, totalLength);
	/* split kmer in several parts */
	}else{
		part = 0;
		fillProbabilityHash = true;
		gapOffset = 0;
		T = 0; /* TODO B&B threshold */
		prefixList.clear();
		set_kmer_to_Pvalue_hash(1, MAX_KMER_SIZE);
		scoreMapArray.clear();

		part = 1;
		fillProbabilityHash = false;

		recalculatePrefixProbabilities(prefixList);

		const int start = MAX_KMER_SIZE+1;
		const int stop = totalLength;
		/* for short kmers use less stringent Threshold */
		T = std::min(0, (int)((-MAX_KMER_SIZE+stop-start+1)*1e6)); /* TODO B&B threshold */
		gapOffset = MAX_KMER_SIZE;
		set_kmer_to_Pvalue_hash(start, stop);
		scoreMapArray.clear();
	}
}

/* calculate probability of prefixes given the PWM starting at the given match column */
template<class HASH>
void Kmer_Generator_Reloaded<HASH>::recalculatePrefixProbabilities(
    prefix_list_t &kmers) {
  float max = 0;
  for (prefix_list_t::iterator it = kmers.begin(); it != kmers.end(); ++it) {
    int col = 1;
    it->probability = 1;
    for (int i = 0; i < it->kmer.length(); ++i) {
      const uint8_t chr = it->kmer.charAt(i);
      if (chr != 0) {
        it->probability = static_cast<float> (it->probability
            * PWM_lin[col][chr]);
        ++col;
      }
    }
    if (it->probability > max) {
      max = it->probability;
    }
  }

  float normalizeProbs = 0;
  for (prefix_list_t::iterator it = kmers.begin(); it != kmers.end();) {
    /* Threshold 0.8 has been chosen as trade-off between performance
     * on benchmarks (Harbison) and E-Value distributions on randomly
     * sampled sequences.
     */
    if (it->probability < 0.8 * max) {
      it = kmers.erase(it);
    } else {
      normalizeProbs += it->probability;
      ++it;
    }
  }

  for (prefix_list_t::iterator it = kmers.begin(); it != kmers.end(); ++it) {
    it->probability /= normalizeProbs;
  }

  rightKmerCond.clear();
  rightKmerCondComp.clear();
}

template<class HASH>
float Kmer_Generator_Reloaded<HASH>::getRightKmerConditional(const UngappedKmer& k){
	if(rightKmerCond.find(k) != rightKmerCond.end()){
		return rightKmerCond[k];
	}
	double condProb = 0;
	int start = MAX_KMER_SIZE - NullModel::getOrder() + k.length() -1;
	int stop = MAX_KMER_SIZE - 1;
// 	cerr << "start: " << start << "\tstop: " << stop << endl;
	for (prefix_list_t::iterator it = prefixList.begin(); it != prefixList.end(); ++it) {
// 			cerr << endl << "suffix: " << k << endl;
// 			cerr << "part0: " << it->kmer << "\t prob: " << it->probability << endl;
// 			cerr << "prefix: " << it->kmer.subKmer(start, stop) << "\t total: " << it->kmer.subKmer(start, stop) * k << endl;
	//		cerr << "conditional: " << NullModel::getConditionals()[it->kmer.subKmer(start, stop) * k] << endl;
		condProb = std::max(condProb, NullModel::getConditional(it->kmer.getLastPosCombined(stop-start+1,k)));
		//condProb += it->probability * NullModel::getConditionals()[it->kmer.subKmer(start, stop) * k];
	}

	rightKmerCond[k] = (float)condProb;
	return (float)condProb;
}



template<class HASH>
Kmer_Generator_Reloaded<HASH>::~Kmer_Generator_Reloaded(){
	if (PWM!=NULL) {
		for (int i = 0; i <= PWM_LENGTH; ++i) {
			delete[] PWM[i];
		}
		delete[] PWM;
	}
}



template<class HASH> template <class Alloc>
inline unsigned int Kmer_Generator_Reloaded<HASH>::binarySearch(const vector<score_vector, Alloc>& a, int S){
	 int low = 0;
	 int high = static_cast<int>(a.size());
	 while (low < high) {
		int mid = low + ((high - low) / 2);
		if (a[mid].score > S)
			low = mid + 1;
		else
			//can't be high = mid-1: here A[mid] >= value,
			//so high can't be < mid if A[mid] == value
			high = mid;
		}
	// high == low, using high or low depends on taste
	if ((low < static_cast<int>(a.size())) && (a[low].score == S))
		return low; // found
	else
		return low-1; // next better score in list
}



template<class HASH>
inline double Kmer_Generator_Reloaded<HASH>::get_kmer_probability_all_positions(unsigned char* kmer) {
	static unsigned char stripped_kmer[100];
	int matchNum = 0;
	int pos = 0;
	while (matchNum < totalLength) {
		stripped_kmer[matchNum] = kmer[pos];
		pos += 1 + gaps[matchNum];
		++matchNum;
	}
	double prob = get_kmer_probability_matched_only(stripped_kmer);
	if (debug) {
		kmer[pos] = '\0';
		stripped_kmer[matchNum] = '\0';
		fprintf(stderr, "g_k_p(%s;%s) = %e\n", translateSeq(kmer,static_cast<int>(strlen((char*)kmer))).c_str(),
				translateSeq(stripped_kmer,static_cast<int>(strlen((char*)stripped_kmer))).c_str(), prob);
	}
	return prob;
}


template<class HASH>
inline double Kmer_Generator_Reloaded<HASH>::get_kmer_probability_matched_only(unsigned char* kmer){
	double pValue = 0;
	id_type id = 0;
	if(totalLength <= MAX_KMER_SIZE){
		for(int i=totalLength-1; i>=0; i--){
			if(kmer[i] == 0) return 1.0;
			id = (id << ID_SHIFT) + kmer[i] - 1;
		}
		//if (debug) std::cerr << "id: " << id << std::endl;
		if(pValueHash[0][id] == 0) return 1.0;
		else return pValueHash[0][id];
	}else{
		/* calculate pVal for first part, return infinity if not in hash */
		for(int i=MAX_KMER_SIZE-1; i>=0; i--){
			if(kmer[i] == 0) return 1.0;
			id = (id << ID_SHIFT) + kmer[i] - 1;
		}
		//cerr << "firstPart: " << pValueHash[0][id] << endl;
		if(pValueHash[0][id] == 0) return 1.0;
		double pVal_0 = pValueHash[0][id];

		/* calculate pVal for second part, return infinity if not in hash */
		id=0;
		for(int i=totalLength-1; i>=MAX_KMER_SIZE; i--){
			if(kmer[i] == 0) return 1.0;
			id = (id << ID_SHIFT) + kmer[i] - 1;
		}
		//cerr << "lastPart: " << pValueHash[1][id] << endl;
		if(pValueHash[1][id] == 0) return 1.0;
		double pVal_1 = pValueHash[1][id];
		//cerr << "\npVal_0 * pVal_1 = " << pVal_0 << "*" << pVal_1 << " = " << pVal_0 * pVal_1 
		//	  << " > " << posSetSize / (totalPositions*2.0 ) << endl;
		//cerr << "posSetSize: " << posSetSize << "\ttotalPositions: " << totalPositions << endl;
		/* if pValues are too bad, return infinity */
		if(pVal_0 * pVal_1 > posSetSize / (totalPositions*2.0 )) return 1.0;

		const int s = get_Kmer_score(kmer);
		const int s1_max = scoreIdVector[0][0].score;
		const int s2_max = scoreIdVector[1][0].score;
		const int vector0_size = static_cast<int>(scoreIdVector[0].size());
		//unsigned int vector1_size = scoreIdVector[1].size();

		/* find next better value in scoreIdVector[1] given score s-s1_max */
		int it_s_s1 = binarySearch(scoreIdVector[1], s-s1_max);

		/* all s1 scores with s1 >= s-s2_max */
		for(int it_s1 = 0 ; it_s1 < vector0_size; it_s1++){
			const int s1 = scoreIdVector[0][it_s1].score;
			if(s-s2_max > s1) break;

			//while(it_s_s1 < vector1_size && scoreIdVector[1][it_s_s1].score < s-s1){
			while(it_s_s1 != 0 && scoreIdVector[1][it_s_s1].score < s-s1){
				it_s_s1--;
			}
			const double pVal_s_s1_max = pValueHash[1][scoreIdVector[1][it_s_s1].id_list.front()];

			const id_list_type& id_vector = scoreIdVector[0][it_s1].id_list;
			

			double leftPartProb = 0;
			for(id_list_type::const_iterator it_vec = id_vector.begin(); it_vec != id_vector.end(); it_vec++){
				leftPartProb += probabilityHash[*it_vec];
			}
			//cerr << "pValue += " << leftPartProb << " * " << pVal_s_s1_max;
			//cerr << " = " << leftPartProb * pVal_s_s1_max << endl;
			pValue += leftPartProb * pVal_s_s1_max;
		}
	}
// 	for (int i=0; i<totalLength; ++i) {
// 		fprintf(stderr, "%c", AlphaChar(kmer[i], Global::A));
// 	}
// 	fprintf(stderr, ",  pValue: %e\n\n", pValue);

	return pValue;
}




/* create a list of k-mers which have a better p-value than a threshold for motif M */
template<class HASH>
void Kmer_Generator_Reloaded<HASH>::set_kmer_to_Pvalue_hash(int start, int stop){
	counts = 0;
	length = stop-start+1;

	unsigned char* idx = (unsigned char*)calloc(length+1, sizeof(unsigned char));

	int i, j, swapIndex;
	/* create the array rowSwap
	 * gives a mapping for the order in which the values appear in every colomn in pwm
	 * eg for first colomn: pwm (3 2 8 5) => rowSwap (3 4 1 2)
	 */
	rowSwap = (int**)calloc(length+1, sizeof(int*));
	for(int i=0; i<=length; i++)rowSwap[i] = (int*)calloc(nAlpha(Global::A)+1, sizeof(int));
	for(i=start; i<=stop; i++){
		for(j=1; j<=nAlpha(Global::A);j++){
			if(j==1) rowSwap[i-start+1][j] = 1;
			else{
				swapIndex = j;
				while(swapIndex > 1 &&
					PWM[i][rowSwap[i-start+1][swapIndex-1]] < PWM[i][j]){
					rowSwap[i-start+1][swapIndex] = rowSwap[i-start+1][swapIndex-1];
					swapIndex--;
				}
				rowSwap[i-start+1][swapIndex] = j;
			}
		}
	}

	/* create the array colSwap
     * gives a mapping for the order in which the values of the highest value of every colomn appear in pwm
	 * eg best value in every row: pwm (3 2 8 5) => colSwap (3 4 1 2)
	 */
	colSwap = (int*)calloc(length+1, sizeof(int));
	for(i=start; i<=stop; i++){
		if(i==start)colSwap[i-start+1] = 1;
		else{
			swapIndex = i-start+1;
			while(swapIndex > 1 &&
				PWM[colSwap[swapIndex-1]][rowSwap[colSwap[swapIndex-1]][1]]<
				PWM[i][rowSwap[i-start+1][1]]){
				colSwap[swapIndex] = colSwap[swapIndex-1];
				swapIndex--;
			}
			colSwap[swapIndex] = i-start+1;
		}
	}

	//fprintf(stderr,"Matrix of %s part\n", part==0 ? "left" : "right");
	//printSwappedPWM(PWM, start, stop);
	//printSwapPWM(rowSwap);
	//fprintf(stderr, "colSwap: \n"); for(i=1; i<=length; i++)fprintf(stderr, "%d\t", colSwap[i]); fprintf(stderr, "\n\n");

	/* create the array swapped PWM with log odds
	 * rearranged array of pwm, in which in every colomn the values are sorted and the
	 * rows are sorted with the highest value in every colomn	 *
	 */
	swappedPWM = (double**)calloc(length+1, sizeof(double*));
	swappedPWM[0] = (double*)calloc(nAlpha(Global::A)+1, sizeof(double));
	for(i=start; i<=stop; i++){
		swappedPWM[i-start+1] = (double*)calloc(nAlpha(Global::A)+1, sizeof(double));
		for(j=1; j<=nAlpha(Global::A); j++){
			swappedPWM[i-start+1][j] = PWM[colSwap[i-start+1]+start-1][rowSwap[colSwap[i-start+1]][j]];
		}
	}

	/* create array S_max
	 * stores the maximally attainable score for positions j to k
	 */
	S_max = (double*)calloc(length+2, sizeof(double));
	for(i=length; i>0; i--){
		if(i < length)S_max[i] = S_max[i+1] + swappedPWM[i][1];
		else S_max[i] = swappedPWM[i][1];
	}

	create_similar_kmers_rec(1,0,idx);

	/* create hash which maps to every kmer (key) the p-value (value) */
	double pValue = 0, pValue_single = 0;


	//cerr << "!!! calculate Pvalues !!!" << endl;
	for(mapType::const_reverse_iterator it = scoreMapArray.rbegin(); it != scoreMapArray.rend(); ++it){
		/* copy hash values to vector */
		scoreIdVector[part].push_back(score_vector(it->first, it->second));
		/* calculate the sum of all pValues for the kmers with the same score */
		const id_list_type& valueList = it->second;
		for(id_list_type::const_iterator it2 = valueList.begin(); it2 != valueList.end(); ++it2){
			//fprintf(stderr, "%lu\n", *it2);
			pValue_single = calculate_kmer_probability(*it2);
			//fprintf(stderr, "pValue: %f\n", pValue_single);
			if(fillProbabilityHash) probabilityHash[*it2] = pValue_single;
			pValue += pValue_single;
		}

		/* insert kmers with pValue into the hash */
		//fprintf(stderr, "%d : ",  it->first);
		for(id_list_type::const_iterator it2 = valueList.begin(); it2 != valueList.end(); ++it2){
			//std::cerr << "\t" << *it2 << ":" << calculate_kmer_probability(*it2) << "=>" << pValue << "\n";
			pValueHash[part][*it2] = pValue;
		}
		//cout << "\n";
	}
	//fprintf(stderr, "threshold: %f, %d different kmers created\n", T*1e-6, counts);
	//if(part == 1)exit(-1);

	/* free allocated memory */
	for(i=0; i<=length; i++)free(rowSwap[i]);
	free(rowSwap);
	for(i=0; i<=length; i++)free(swappedPWM[i]);
	free(swappedPWM);
	free(colSwap);
	free(idx);
	free(S_max);
}



/* j: current position in the generated k-mer
 * S_j: score for the generated k-mer up to position j-1
 * idx: string for the found k-mer
 */
template<class HASH>
void Kmer_Generator_Reloaded<HASH>::create_similar_kmers_rec(int j, double S_j, unsigned char* idx){
	double S_j_1;
	for(int a=1;a<=nAlpha(Global::A);a++){
		S_j_1 = S_j + swappedPWM[j][a];
		if((int)((S_j_1 + S_max[j+1])*1e6) < T) continue;
		idx[colSwap[j]-1] = static_cast<unsigned char>(rowSwap[colSwap[j]][a]);  // store original nucleotide into string
		if(j+1 <= length) create_similar_kmers_rec(j+1, S_j_1, idx);
		else{
			idx[j] = '\0';

			int intScore = (int)(S_j_1 * 1e6); // use integer value as key to avoid rounding errors due to double precision
			id_type id =0;
			for(int pos=length-1; pos>=0; pos--) id = (id << ID_SHIFT) + idx[pos] - 1;
			scoreMapArray[intScore].push_back(id); // include string into list
			counts++;
		}
	}
}



template<class HASH>
void Kmer_Generator_Reloaded<HASH>::printSwapPWM(int **pwm){
	int i, j;
	fprintf(stderr, "\n\t");
		for(i=0;i<length;i++) fprintf(stderr, "%6d\t", i+1);
		fprintf(stderr, "\n--------");
		for(i=0;i<length;i++)fprintf(stderr, "--------");
		for(i=0;i<=nAlpha(Global::A);i++){
		 	if(i==0){
		 		fprintf(stderr, "\ncounts\t");
		 		for(j=1; j<=length; j++){
		 			fprintf(stderr,"%d\t", pwm[j][i]);
		 		}
		 		fprintf(stderr, "\n");
		 		continue;
		 	}
		 	fprintf(stderr, "\n%c\t", AlphaChar(i, Global::A));
		 	for(j=1; j<=length; j++){
		 		fprintf(stderr,"%d\t", pwm[j][i]);
		 	}
		}
		fprintf(stderr, "\n\n");
}



template<class HASH>
void Kmer_Generator_Reloaded<HASH>::printSwappedPWM(double **pwm, const int start, const int stop){
	int i, j;
	fprintf(stderr, "\n\t");
		for(i=0;i<totalLength;i++) {
			if (i!=0 && (i<start || i>stop)) continue;
			fprintf(stderr, "%6d\t", i+1);
		}
		fprintf(stderr, "\n--------");
		for(i=0;i<totalLength;i++) {
			if (i!=0 && (i<start || i>stop)) continue;
			fprintf(stderr, "--------");
		}
		for(i=0;i<=nAlpha(Global::A);i++){
		 	if(i==0){
		 		fprintf(stderr, "\ncounts\t");
		 		for(j=1; j<=totalLength; j++){
					if (j<start || j>stop) continue;
		 			fprintf(stderr,"%6.0f\t", pwm[j][i]);
		 		}
		 		fprintf(stderr, "\n");
		 		continue;
		 	}
		 	fprintf(stderr, "\n%c\t", AlphaChar(i, Global::A));
		 	for(j=1; j<=totalLength; j++){
				if (j<start || j>stop) continue;
		 		fprintf(stderr,"%6.2f\t", pwm[j][i]);
		 	}
		}
		fprintf(stderr, "\n\n");
}

#endif
