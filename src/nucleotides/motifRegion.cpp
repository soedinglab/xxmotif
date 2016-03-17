#include "motifRegion.h"

MotifRegion::MotifRegion(int startRange, int endRange, double enrichmentThreshold, double significanceThreshold) :
		_startRange(startRange),
		_endRange(endRange),
		_enrichmentThreshold(enrichmentThreshold),
		_significanceThreshold(significanceThreshold)
		{
	LOG_Bonferonni = LogTable::LOG_i[endRange] + LogTable::LOG_i[endRange-1] - LogTable::LOG_i[2];
}

/* calculate max of distribution by using a sliding window of size windowSize */
int MotifRegion::findMaximum(int* startPosCounts){
	int max = 0;
	double maxOccurrence = 0;	// highest sliding window score

	// initialize positions array with -1
	for(int i=0; i<=_windowSize; i++)_positions[i] = -1;

	const int halfWS = _windowSize / 2;

	for(int j = _startRange; j< _endRange; j++){
		/* store positions if not all previous positions are known yet */
		if(_positions[_windowSize-1] == -1){
			int pos = 0;
			/* find first unknown position */
			while(_positions[pos] != -1) pos++;
			_positions[pos] = j;
		}else{
			_positions[_windowSize] = j;
			double occ = 0.0;
			int current = _positions[halfWS];
			/* calculate window score */
			for(int i=0; i<halfWS; i++){
				occ += startPosCounts[_positions[i]]/(current - _positions[i] + 1);
			}
			occ += startPosCounts[current];
			for(int i=halfWS + 1; i<= _windowSize; i++){
				occ += startPosCounts[_positions[i]]/(_positions[i] - current + 1);
			}

			if(occ > maxOccurrence){
				max = current;
				maxOccurrence = occ;
			}
			/* shift window */
			for(int i=0;i<_windowSize;i++){
				_positions[i] = _positions[i+1];
			}
		}
	}

	/* if less than windoSize position with motif, use the position with the most occurrences as max */
    if(_positions[_windowSize] == -1){
	    int maxNb = 0;
	    for(int i=0; i<_windowSize; i++){
	    	if(_positions[i] == -1)	break;
	    	else if(startPosCounts[_positions[i]] > maxNb){
	    	    max = _positions[i];
	    	    maxNb = startPosCounts[_positions[i]];
	    	}
	    }
	}
    return max;
}

region MotifRegion::getRegion(int* startPosCounts, int counts){
	int max = findMaximum(startPosCounts); 			// position with highest sliding window score
	int L = _endRange-_startRange+1; // positions with possible motif

	//fprintf(stderr, "\n\nmax: %d, possible motif positions: %d\n", max, L);
	double newProb;
	double oldBestProb = 0;

	/* initialize enriched region to only the maximum */
	int x1 = max; 		// left end of region
	int x2 = max; 		// right end of region

	/* calculate number of motifs in enriched region */
	int nRound = 0;
	for(int i=x1; i<=x2; i++) nRound += startPosCounts[i];

	int occ = nRound;

	/* calculate probability of motif overrepresentation by normal approximation of binomial distribution */
	double bestProb = calcRegionPval(x1, x2, L, occ, counts);

	if(occ*1.0/counts > _enrichmentThreshold){
		oldBestProb = bestProb;
	}else{
		oldBestProb = 100;
	}

	//fprintf(stderr, "x1: %d, x2: %d, log z: %f, z: %e, occ: %d, counts: %d\n", x1, x2, bestProb, exp(bestProb), occ, counts);

	/* extend region to both sides */
	int x1Prime = x1;
	int loop, x1Test, x2Test;
	do{
		loop = 0;

		bestProb = oldBestProb;
		/* extend region to the left */
		for(x1Test = x1-1; x1Test >= _startRange; x1Test--){
			if(startPosCounts[x1Test] == 0) continue;
			occ += startPosCounts[x1Test];
			newProb = calcRegionPval(x1Test, x2, L, occ, counts);

			//fprintf(stderr, "x1: %d, x2: %d, log z: %f, z: %e, occ: %d, counts: %d\n", x1Test, x2, newProb, exp(newProb), occ, counts);
			if(newProb < bestProb){
				bestProb = newProb;
				x1Prime = x1Test;
			}
		}
		double bestProbLeft = bestProb;
		//fprintf(stderr, "bestProbLeft: %e => %e, x1Prime: %d\n", exp(bestProbLeft), exp(bestProbLeft+LOG_Bonferonni), x1Prime);
		int x2Prime = x2;

		/* reset bestProb and motif occurrencies to the beginning of the round */
		occ = nRound;
		bestProb = oldBestProb;

		/* extend region to the right */
		for(x2Test = x2+1; x2Test <= _endRange; x2Test++){
			if(startPosCounts[x2Test] == 0) continue;
			occ += startPosCounts[x2Test];
			newProb = calcRegionPval(x1, x2Test, L, occ, counts);

			//fprintf(stderr, "x1: %d, x2: %d,log z: %f, z: %e, occ: %f, counts: %d\n", x1, x2Test, newProb, exp(newProb), occ, counts);
			if(newProb < bestProb){
				bestProb = newProb;
				x2Prime = x2Test;
			}
		}
		occ = nRound;
		double bestProbRight = bestProb;

		LOG(DEBUG4) << strprintf("bestProbRight: %e => %e, x2Prime: %d", exp(bestProbRight), exp(bestProbRight+LOG_Bonferonni),  x2Prime);

		LOG(DEBUG4) << strprintf("\nbestProbLeft: %e, bestProbRight: %e, oldBestProb: %e, x1: %d, x1Prime: %d, x2: %d, x2Prime: %d", exp(bestProbLeft), exp(bestProbRight), exp(oldBestProb), x1,x1Prime,x2, x2Prime);

		if(std::min(bestProbLeft, bestProbRight) < oldBestProb){
			if(oldBestProb == 100){
				if(bestProbLeft < bestProbRight){
					x1 = x1Prime;
				}else{
					x2 = x2Prime;
				}
			}else{
				/* combine best region ends of both sides */
				x1 = x1Prime;  // best region end left
				x2 = x2Prime;  // best region end right
			}
			nRound = 0;
			for(int i=x1;i<=x2;i++)nRound += startPosCounts[i];
			occ = nRound;
			/* only accept region when it has at least 20 % of all occurrences */
			bestProb = calcRegionPval(x1, x2, L, occ, counts);
			if(occ*1.0/counts > _enrichmentThreshold){
				oldBestProb = bestProb;
			}else{
				oldBestProb = 100;
			}
			LOG(DEBUG4) << strprintf("\noldBestProb: x1: %d, x2: %d, log z: %f, z: %e, occ: %f, counts: %d", x1, x2, oldBestProb, exp(oldBestProb), occ, counts);

			loop = 1;
		}
	}while(loop);
	/* possible start end situations: seqLength over 2 */

	bestProb += LOG_Bonferonni;

	region r;
	if(occ*1.0/counts > _enrichmentThreshold && bestProb < _significanceThreshold){
		LOG(DEBUG4) << strprintf("\n\nenrichment: %f, occ: %f, counts: %d, bestProb: %f, x1: %d, x2: %d\n", occ*1.0/counts, occ, counts, bestProb, x1, x2);
		r.set = 1;
		r.max = max;
		r.startRegion = x1;
		r.endRegion = x2;
	}else{
		LOG(DEBUG4) << strprintf("\n\nrange1: %d, range2: %d\n", _startRange, _endRange);
		r.set = 0;
		r.max = _startRange;
		r.startRegion = _startRange;
		r.endRegion = _endRange;
	}
	return r;
}
