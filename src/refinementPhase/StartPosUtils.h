#ifndef STARTPOS_UTILS
#define STARTPOS_UTILS

#include "Sorted_Sites.h"
#include "../Globals.h"

bool updateMultOcc(sorted_sites* sites, int& sites_counter, int32_t seq, int32_t pos, int length, float pVal, ss_type set);

void sortSitesByScore(sorted_sites* sites, int nbSites);
void sortSitesByStartPos(sorted_sites* sites, int nbSites, bool palin=false);

inline bool updateMultOcc(sorted_sites* sites, int& sites_counter, int32_t seq, int32_t pos, int length, double pVal, ss_type set){

	/* detect overlapping motifs */
	if(sites_counter > 0 &&
	  sites[sites_counter-1]->sequence == seq &&
	  sites[sites_counter-1]->startPos+length >= pos){
		/* if last match has better pVal, keep it */
		if(sites[sites_counter-1]->score <= pVal){
			return false;
		}
		/* if new match has better pVal, overwrite last match */
		else sites_counter--;
	}


	//if(debug)fprintf(stderr, "seq: %d, pos: %d, pVal: %e, pValScore: %e, bestSeqNb: %d, pValCons: %e\n\n", seq, pos, pVal, pValScore, bestSeqNb, pValCons);
	/* remove overlapping motifs on reverse strand */
	if(Global::revcomp && pos > length / 2){
		int back=1;
		int revPos = set->entity[seq]->n - length + 2 - pos;
		while(sites_counter - back >= 0 &&
			  sites[sites_counter - back]->sequence == seq &&
			  sites[sites_counter - back]->startPos > revPos + length -1){
			back++;
		}
		/* if there is an overlap on the reverse strand */
		if(sites_counter - back >= 0 &&
		   sites[sites_counter - back]->sequence == seq &&
		   sites[sites_counter - back]->startPos >= revPos -length +1){

				/* if previous hit is better then current hit , don't use current hit*/
			if(sites[sites_counter - back]->score <= pVal){
				return false;
			}
			/* if new match has better pVal, overwrite last match */
			else{
				while(back > 1){
					sites[sites_counter - back]->sequence = sites[sites_counter - back + 1]->sequence;
					sites[sites_counter - back]->startPos = sites[sites_counter - back + 1]->startPos;
					sites[sites_counter - back]->score = sites[sites_counter - back + 1]->score;

					back--;
				}
				sites_counter--;
			}
		}
	}

	sites[sites_counter]->sequence = seq;
	sites[sites_counter]->startPos = pos;
	sites[sites_counter]->score = static_cast<float>(pVal);

	return true;
}


inline void sortSitesByScore(sorted_sites* sites, int nbSites){
	qsort (sites, nbSites, sizeof (sorted_sites), SortedSites::compare_scores);
}

inline void sortSitesByStartPos(sorted_sites* sites, int nbSites, bool palin){
	if(palin){
		qsort (sites, nbSites, sizeof (sorted_sites), SortedSites::compare_startPos_palin);
	}else{
		qsort (sites, nbSites, sizeof (sorted_sites), SortedSites::compare_startPos);
	}
}

#endif
