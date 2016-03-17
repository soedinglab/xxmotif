#ifndef ELONG_CANDIDATES_H_
#define ELONG_CANDIDATES_H_

#include <unordered_set>
#include <memory>

#include "startPosList.h"
#include "Kmer.h"
#include "Extension.h"
#include "../SmallKmer.h"
#include "../ProgressBar.h"
#include "../Globals.h"

/* class providing all interesting kmers which should be elongated in the next step */
class ElongCandidates {

public:
	ElongCandidates(int gapIndex, motif_type motType, double pValThreshold);
	~ElongCandidates();
	void getElongationCandidates(elongList& kmerList);

	static bool compareKmerResults(std::shared_ptr<Kmer>& a, std::shared_ptr<Kmer>& b);
	static void freeMemory();
	//static bool compareKmerResults(Kmer& a, Kmer& b);

	static const char IUPAC_CHARS[15];
	static const uint8_t IUPAC_palindrome[11];
	static double LOG_Bonferonni;
	static int nbSequences;
	static int totalPositions;
	static double avgLength;

	static int* motifsPerSequenceCount;	/* stores how many motifs are found per sequence */
	static motif_type type;

private:
	static fullHashType *extendedMotifsHash;
	static uint32_t *countIUPAC;		/* number of occurrences of kmer */
	static double   *kmerProbs;			/* pValue for every kmer */

	static bool initialized;

	startPosList* collectStartPosHash;
	uint8_t _gapIndex;

	void countAllOccurrences();

	template <typename elemType>
	void mergeIUPACS_rec(elemType* countIUPAC, uint32_t id, elemType& value, int pos, int depth);

	void createIUPAC_list_phase1_rec(std::list<std::pair<uint32_t, double> >& bestPhaseList, int& listSize,
			uint32_t minOcc, int N, uint32_t id, int pos, int depth);

	void createIUPAC_list_phase2_rec(std::list<std::pair<uint32_t, double> >& bestPhaseList, int& listSize,
				int minOcc, int N, uint32_t id, MotifRegion& motifRegion, int pos, int depth);

	void setKmerProbs();
	int setLog_Pvalues();

	void createStartPositionList();
	void setStartPos();

	int filterElongationCandidates(elongList& kmerList);
	double getLog_Pvalue(uint32_t id, startPosList& spl, int startRange, int endRange, bool mops);

	uint32_t StringToId(const char* s);
	void printString(uint32_t id);

	void updateProgressBar(int step);
	void setTrackedMotif(const std::unordered_set<std::string>& motif);


	/* ********************************
	 * CONSTANTS
	 * ********************************/
	static int eraseBits[];						/* bit mask for erasing bits in uint32_t */

	static uint32_t POW_16_coding;
	static const uint32_t POW_16_3 = 4096;

	static uint32_t POW_4_coding;
	static const uint32_t POW_4_3 = 64;

	static const uint8_t max_pattern_width = 22;    /* size of kmer in which the coding characters have to be located */
	static uint8_t coding_characters;   			/* number of coding characters in kmer */
	static const uint8_t IUPAC_characters = 2;      /* maximal number of iupac characters in kmer */
	static const uint8_t IUPAC_CHARS_NB = 10;	    /* states: A,C,G,T,M,R,W,S,Y,K */

	static const unsigned int NUMBER_PHASE1_DEGENERATIONS;
	static const unsigned int NUMBER_PHASE2_DEGENERATIONS;
	static const unsigned int NUMBER_PHASE3_DEGENERATIONS;

	static const uint8_t IUPACS[4][5];

#ifdef SIXMER
	static const uint8_t gaps[56][5];
#else
	static const uint8_t gaps[35][4];
#endif

	static const uint8_t gaps_palindrome[16][5];

	static std::list< std::pair<uint8_t,uint32_t> > TEST_ID;
};

template <typename elemType >
void ElongCandidates::mergeIUPACS_rec(elemType* array, uint32_t id, elemType& value, int pos, int depth){
	uint32_t newid, tmp = id;
	/* increment number of IUPAC characters in current string */
	depth++;
	/* shift to pos position */
	tmp = tmp >> pos*4;

	/* every position which is to degenerate */
	for (int i = pos; i < coding_characters; i++) {
		const uint8_t *IUPAC_base = IUPACS[tmp & 15];
		tmp = tmp >> 4;

		/* null out ith character */
		uint32_t erased_id = id & eraseBits[i];

		/* every allowed degeneration */
		for(int j=1; j<4; j++){
			/* add new iupac character at ith position */
			newid = erased_id | (IUPAC_base[j] << 4*i);

			/* combine information */
			array[newid] += value;

			if(depth < IUPAC_characters){
				/* recursivly go to more degenerated IUPAC strings */
				mergeIUPACS_rec(array, newid, value, i+1, depth);
			}
		}
	}
}

inline bool ElongCandidates::compareKmerResults(std::shared_ptr<Kmer>& a, std::shared_ptr<Kmer>& b){
	if(a->p_set < b->p_set) return true;
	else return false;
}

inline void ElongCandidates::freeMemory(){
	free(countIUPAC);
	free(kmerProbs);
	free(motifsPerSequenceCount);
	delete extendedMotifsHash;
}

inline void ElongCandidates::updateProgressBar(int step){
	if(!Global::batch_mode){
		if(type == FIVEMERS){
#ifdef SIXMER
			cout << "\relongate sixmers:              "<< progressBar(6*_gapIndex+step, 6*Global::GAPS-1, 40) << std::flush;
#else
			cout << "\relongate fivemers:             "<< progressBar(6*_gapIndex+step, 6*Global::GAPS-1, 40) << std::flush;
#endif
		}else if(type == PALINDROME){
			cout << "\relongate palindromic sixmers:  "<< progressBar(6*_gapIndex+step, 6*(Global::maxMatchPositions-5)-1, 40) << std::flush;
		}else if(type == TANDEM){
			cout << "\relongate tandemic sixmers:     "<< progressBar(6*_gapIndex+step, 6*(Global::maxMatchPositions-5)-1, 40) << std::flush;
		}
	}
}

#endif /* ELONG_CANDIDATES_H_ */
