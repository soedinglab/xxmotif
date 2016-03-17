#include "elongationCandidates.h"
#include "../NullModel.h"
#include "../utils.h"
#include "../nucleotides/motifRegion.h"
#include "../AbstractKmer.h"
#include "../Globals.h"
#include "../SmallKmer.h"
#include "../elongationPhase/Kmer.h"
#include "../elongationPhase/Match.h"
#include "../elongationPhase/elongationCandidates.h"
#include "../memoryPool/pool_alloc.h"

#include <iomanip>
#include <fstream>

bool ElongCandidates::initialized = false;
int ElongCandidates::eraseBits[6];
int ElongCandidates::nbSequences;
int ElongCandidates::totalPositions;
double ElongCandidates::avgLength;

fullHashType *ElongCandidates::extendedMotifsHash;
uint32_t *ElongCandidates::countIUPAC;			/* number of occurrences of kmer */
double   *ElongCandidates::kmerProbs;			/* pValue for every kmer */
int		 *ElongCandidates::motifsPerSequenceCount;

motif_type ElongCandidates::type = FIVEMERS;

/*************
 * CONSTANTS
 *************/
/* number of degenerations per kmer that should be passed to localization pvalue calculation*/
const unsigned int ElongCandidates::NUMBER_PHASE1_DEGENERATIONS = 20;
/* number of degenerations per kmer that should be passed to elongation*/
const unsigned int ElongCandidates::NUMBER_PHASE2_DEGENERATIONS = 5;
/* number of extended degenerations per kmer that should be passed to pwm refinement*/
const unsigned int ElongCandidates::NUMBER_PHASE3_DEGENERATIONS = 5;
uint8_t ElongCandidates::coding_characters;
uint32_t ElongCandidates::POW_16_coding;
uint32_t ElongCandidates::POW_4_coding;
double ElongCandidates::LOG_Bonferonni;


const uint8_t ElongCandidates::IUPACS[4][5] = { { 0, 4, 5, 6, 10 }, /* A, M, R, W, N */
							          { 1, 4, 7, 8, 10 }, /* C, M, S, Y, N */
							          { 2, 5, 7, 9, 10 }, /* G, R, S, K, N */
							          { 3, 6, 8, 9, 10 } }; /* T, W, Y, K, N */

												// A C G T M R W S Y K N
const uint8_t ElongCandidates::IUPAC_palindrome[11] = {3,2,1,0,9,8,6,7,5,4,10};
											  // 0    1    2    3    4    5    6    7    8    9    10   11   12   13   14
const char ElongCandidates::IUPAC_CHARS[15] = { 'A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'N', 'B', 'D', 'H', 'V' };

#ifdef SIXMER
const uint8_t ElongCandidates::gaps[56][5] = {
	/* 0 gaps */ {0,0,0,0,0},

	/* 1 gap  */ {1,1,1,1,1},{0,1,1,1,1},{0,0,1,1,1},{0,0,0,1,1},{0,0,0,0,1},

	/* 2 gaps */ {2,2,2,2,2},{1,2,2,2,2},{1,1,2,2,2},{1,1,1,2,2},{1,1,1,1,2},
				 {0,2,2,2,2},{0,1,2,2,2},{0,1,1,2,2},{0,1,1,1,2},{0,0,2,2,2},
				 {0,0,1,2,2},{0,0,1,1,2},{0,0,0,2,2},{0,0,0,1,2},{0,0,0,0,2},

	/* 3 gaps */ {3,3,3,3,3},{2,3,3,3,3},{2,2,3,3,3},{2,2,2,3,3},{2,2,2,2,3},
				 {1,3,3,3,3},{1,2,3,3,3},{1,2,2,3,3},{1,2,2,2,3},{1,1,3,3,3},
				 {1,1,2,3,3},{1,1,2,2,3},{1,1,1,3,3},{1,1,1,2,3},{1,1,1,1,3},
				 {0,3,3,3,3},{0,2,3,3,3},{0,2,2,3,3},{0,2,2,2,3},{0,1,3,3,3},
				 {0,1,2,3,3},{0,1,2,2,3},{0,1,1,3,3},{0,1,1,2,3},{0,1,1,1,3},
				 {0,0,3,3,3},{0,0,2,3,3},{0,0,2,2,3},{0,0,1,3,3},{0,0,1,2,3},
				 {0,0,1,1,3},{0,0,0,3,3},{0,0,0,2,3},{0,0,0,1,3},{0,0,0,0,3}
	};
#else
const uint8_t ElongCandidates::gaps[35][4] = {
	/* 0 gaps */ {0,0,0,0},

	/* 1 gap  */ {1,1,1,1},{0,1,1,1},{0,0,1,1},{0,0,0,1},

	/* 2 gaps */ {2,2,2,2},{1,2,2,2},{1,1,2,2},{1,1,1,2},
				 {0,2,2,2},{0,1,2,2},{0,1,1,2},{0,0,2,2},
				 {0,0,1,2},{0,0,0,2},

	/* 3 gaps */ {3,3,3,3},{2,3,3,3},{2,2,3,3},{2,2,2,3},
				 {1,3,3,3},{1,2,3,3},{1,2,2,3},{1,1,3,3},
				 {1,1,2,3},{1,1,1,3},{0,3,3,3},{0,2,3,3},
				 {0,2,2,3},{0,1,3,3},{0,1,2,3},{0,1,1,3},
				 {0,0,3,3},{0,0,2,3},{0,0,1,3},{0,0,0,3}
	};
#endif

const uint8_t ElongCandidates::gaps_palindrome[16][5] = {
		{0,0,0,0,0},{0,0,1,1,1},{0,0,2,2,2},{0,0,3,3,3},{0,0,4,4,4},{0,0,5,5,5},
		{0,0,6,6,6},{0,0,7,7,7},{0,0,8,8,8},{0,0,9,9,9},{0,0,10,10,10},{0,0,11,11,11},
		{0,0,12,12,12},{0,0,13,13,13},{0,0,14,14,14},{0,0,15,15,15}
	};

std::list< std::pair<uint8_t,uint32_t> > ElongCandidates::TEST_ID;


ElongCandidates::ElongCandidates(int gapIndex, motif_type motType, double pValThreshold){
	uint8_t pal_chars = 6;
#ifdef SIXMER
	uint8_t non_pal_chars = 6;
#else
	uint8_t non_pal_chars = 5;
#endif
	if(!initialized){
		type = motType;
		if(type == PALINDROME || type == TANDEM){
			coding_characters = pal_chars;
		}else{
			coding_characters = non_pal_chars;
		}
		for(int i=0; i<6; i++){
			eraseBits[i] = 15;
			eraseBits[i] = (eraseBits[i] << (4*i));
			eraseBits[i] = ~eraseBits[i];
		}

		if(!Global::trackedElongation) setTrackedMotif(Global::trackedMotifs);

		nbSequences = Global::posSet->nent;
		totalPositions = Global::posSet->total[nbSequences];
		avgLength = Global::posSet->avgLength[nbSequences];

		int cbits = 4;
		int gbits = 2;
		if(type == PALINDROME || type == TANDEM){
			gbits = 4;
		}
		SmallKmer::init(cbits, gbits);

		POW_16_coding = static_cast<uint32_t>(pow(16, coding_characters));
		POW_4_coding = static_cast<uint32_t>(pow(4, coding_characters));


		extendedMotifsHash = new fullHashType;

		//extendedMotifsHash->rehash(1e7);

		countIUPAC = (uint32_t*)calloc(POW_16_coding, sizeof(uint32_t)); /* number of occurrences of kmer */
		kmerProbs  = (double*)calloc(POW_16_coding, sizeof(double)); /* pValue for every kmer */

		motifsPerSequenceCount = (int*)malloc((Global::posSet->nent + 1)*sizeof(int));

		/* initialize data structures for motif extension */
		Extension::initialize();

		initialized = true;
	}else if(motType != type){
		type = motType;

		if(type == PALINDROME || type == TANDEM){
			coding_characters = 6;
		}else{
			coding_characters = 5;
		}

		if(!Global::trackedElongation) setTrackedMotif(Global::trackedMotifs);

		int cbits = 4;
		int gbits = 2;
		if(type == PALINDROME || type == TANDEM){
			gbits = 4;
		}

		SmallKmer::init(cbits, gbits);

		POW_16_coding = static_cast<uint32_t>(pow(16, coding_characters));
		POW_4_coding = static_cast<uint32_t>(pow(4, coding_characters));

		if(pal_chars != non_pal_chars){
			free(countIUPAC);
			countIUPAC = (uint32_t*)calloc(POW_16_coding, sizeof(uint32_t)); /* number of occurrences of kmer */
			free(kmerProbs);
			kmerProbs  = (double*)calloc(POW_16_coding, sizeof(double)); /* pValue for every kmer */
		}
	}

	collectStartPosHash = new startPosList[POW_16_coding];

	_gapIndex = static_cast<uint8_t>(gapIndex);

	motif_columns_type motif_columns;
	const uint8_t* gap;
	if(type == FIVEMERS) gap = gaps[_gapIndex];
	else gap = gaps_palindrome[_gapIndex];

	int pos = 0;
	motif_columns.push_back(pos++);
	for(int i=1; i<5; i++){
		pos += (gap[i]-gap[i-1]);
		motif_columns.push_back(pos++);
	}

	LOG_Bonferonni = calculate_log_bonferonni(motif_columns, LogTable::LOG_i[Global::neff_discrete]);
}

ElongCandidates::~ElongCandidates() {
	memset(countIUPAC, 0, POW_16_coding*sizeof(uint32_t));
	memset(kmerProbs, 0, POW_16_coding*sizeof(double));

	delete[] collectStartPosHash;

	if(extendedMotifsHash != NULL) extendedMotifsHash->clear();
}


void ElongCandidates::getElongationCandidates(elongList& elongCandidates) {
	updateProgressBar(0);

	/* count how many instances of a kmer with gap kombination k exist */
	countAllOccurrences();
	updateProgressBar(1);

	/* calculate probability for finding any IUPAC string in random sequences with given order */
	setKmerProbs();
	updateProgressBar(2);

	/* calculate probabilities for finding every kmer in at least countIUPAC random sequences */
	setLog_Pvalues();
	updateProgressBar(3);
	//cerr << "count after Phase 1: " << count << endl;

	/* create a list with all start positions of the instances */
	createStartPositionList();
	updateProgressBar(4);

	/* elongate motifs of 5 best degenerative motif for each nondegenerative motif and filter to 1e10  */
	filterElongationCandidates(elongCandidates);
	updateProgressBar(5);
	//cerr << "count after Phase 2: " << count << endl;
}

void ElongCandidates::countAllOccurrences() {
	const uint8_t *gap = gaps[_gapIndex];
	if(type == PALINDROME || type == TANDEM) gap = gaps_palindrome[_gapIndex];

	/* every sequence */
	for (int i = 1; i <= nbSequences; i++) {
		unsigned char *sequence = Global::posSet->entity[i]->S[0];
		/* every startPos */
		//cerr << endl << i << "/" << Global::posSet->entity[i]->n << "/" << Global::posSet->entity[i]->n-coding_characters-gap[coding_characters-2]+1 << endl;
		int seqLen = Global::posSet->entity[i]->n-coding_characters-gap[coding_characters-2]+1;
		if(type == PALINDROME && Global::revcomp) seqLen /= 2;
		for (int j = 1; j <= seqLen; j++) {
			int base = sequence[j]-1;
			if (base < 0) continue;

			bool N_char = false;
			uint32_t id = base;
			for (int l = 1; l < coding_characters; l++) {
				int base = sequence[j+l+gap[l-1]]-1;
				if (base < 0){
					N_char = true;
					break;
				}
				id = (id << 4) + sequence[j+l+gap[l-1]]-1;
			}

			if(!N_char){
				countIUPAC[id]++;
			}
		}
	}

	/* fill IUPAC degenerations with counts */
	for(uint32_t id=0; id < POW_4_coding; id++){
		int tmp = id;
		uint32_t newId = 0;
		for(int j=0; j<coding_characters; j++){
			int base = tmp&3; tmp = tmp >> 2;
			newId = (newId << 4) + base;
		}
		mergeIUPACS_rec(countIUPAC, newId, countIUPAC[newId], 0, 0);
	}
}

/*
 * calculate the probability of a k-mer on the background set with given order
 */

void ElongCandidates::setKmerProbs(){
	static unsigned char expanded_kmer[100];

	const uint8_t *gap = gaps[_gapIndex];
	if(type == PALINDROME || type == TANDEM) gap = gaps_palindrome[_gapIndex];

	int totalPos = coding_characters + gap[coding_characters-2];
	for(uint32_t id=0; id < POW_4_coding; id++){
		/* create gapped kmer for every non degenerated kmer */
		int tmpId = id;
		int pos = totalPos-1;
		uint8_t base = (tmpId & 3);
		expanded_kmer[pos--] = static_cast<uint8_t>(base+1);

		for (int i = 1; i < coding_characters; i++) {
			int gapSize = gap[coding_characters-i-1];
			if(coding_characters-i-2 >= 0) gapSize -= gap[coding_characters-i-2];

			for (int m = 0; m < gapSize; m++) {
				expanded_kmer[pos--] = 0;
			}
			tmpId = (tmpId >> 2);
			base = (tmpId & 3);
			expanded_kmer[pos--] = static_cast<uint8_t>(base+1);
		}

		uint32_t new_id = expanded_kmer[0]-1;
		for(int i=1; i<totalPos; i++){
			if(expanded_kmer[i] == 0)continue;
			new_id = (new_id << 4 ) + expanded_kmer[i]-1;
		}

		double prob = NullModel::getProbability(expanded_kmer, totalPos);
		kmerProbs[new_id] = prob;

		/* add probability to every IUPAC degeneration of the kmer */
		mergeIUPACS_rec(kmerProbs, new_id, prob, 0, 0);
	}
	//printString(TEST_ID.first, TEST_ID.second);
	//fprintf(stderr, "gap: %d, id. %d, kmer Prob: %f\n", gapIndex, TEST_ID.first, kmerProbs[TEST_ID.second]);
}


void ElongCandidates::createIUPAC_list_phase1_rec(std::list<std::pair<uint32_t, double> >& bestPhaseList,
		int& listSize, uint32_t minOcc, int N, uint32_t id, int pos, int depth){
	uint32_t newId, tmp = id;

	/* increment number of IUPAC characters in current string */
	depth++;
	/* shift to pos position */
	tmp = tmp >> pos*4;

	/* every position which is to degenerate */
	int changePositions = coding_characters;
	int IUPAC_chars = IUPAC_characters;
	if(type == PALINDROME || type == TANDEM){
		changePositions = (coding_characters + 1) / 2;
		IUPAC_chars = (IUPAC_characters + 1) / 2;
	}

	for (int i = pos; i < changePositions; i++) {
		const uint8_t *IUPAC_base = IUPACS[tmp & 15];
		tmp = tmp >> 4;

		/* null out ith character */
		uint32_t erased_id = id & eraseBits[i];
		if(type == PALINDROME){
			erased_id = erased_id & eraseBits[coding_characters-i-1];
		}else if(type == TANDEM){
			if(i > 3) erased_id = erased_id & eraseBits[i-3];
			else	  erased_id = erased_id & eraseBits[i+3];
		}
		/* every allowed degeneration */
		for(int j=1; j<4; j++){
			/* add new iupac character at ith position */
			int base = IUPAC_base[j];
			newId = erased_id | (base << 4*i);
			if(type == PALINDROME){
				newId = newId | (IUPAC_palindrome[base] << 4*(coding_characters-i-1));
			}else if(type == TANDEM){
				if(i>3) newId = newId | (base << 4*(i-3));
				else    newId = newId | (base << 4*(i+3));
			}
			if(countIUPAC[newId] > minOcc){
				double p = kmerProbs[newId];

				/* calculate pValue for finding countIUPAC[id] number of instances in set */
				double set_p = calculateOrderStatisticsPvalue(countIUPAC[newId], N, p) + LOG_Bonferonni;

				for(std::list<std::pair<uint8_t,uint32_t> >::iterator it = TEST_ID.begin(); it != TEST_ID.end(); it++){
					//if(newId == StringToId("CACGTG")){
					if(_gapIndex == (*it).first && newId == (*it).second){
						fprintf(stderr, "\n+++++++++++ ");
						printString(newId);
						fprintf(stderr, "\tgap: %d, counts: %d, N: %d, p: %f\t, E-Value: %e\n", _gapIndex, countIUPAC[newId], N, kmerProbs[newId], exp(set_p));
					}
				}

				/* not enough motifs are in best motif list => insert them at correct position */
				if(listSize < (int)NUMBER_PHASE1_DEGENERATIONS){
					std::list<std::pair<uint32_t, double> >::iterator it;
					for(it = bestPhaseList.begin(); it != bestPhaseList.end(); it++){
						if(it->second > set_p)break;
					}
					bestPhaseList.insert(it, std::pair<uint32_t, double>(newId, set_p));
					listSize++;
					/* if motif has better pval than worst motif in best motif list, drop worst element and include new one
					 * at correct position */
				}else if(set_p < bestPhaseList.back().second){
					bestPhaseList.pop_back();
					std::list<std::pair<uint32_t, double> >::iterator it;
					for(it = bestPhaseList.begin(); it != bestPhaseList.end(); it++){
						if(it->second > set_p)break;
					}
					bestPhaseList.insert(it, std::pair<uint32_t, double>(newId, set_p));
				}
			}else{
				for(std::list<std::pair<uint8_t,uint32_t> >::iterator it = TEST_ID.begin(); it != TEST_ID.end(); it++){
					if(_gapIndex == (*it).first && newId == (*it).second){
						fprintf(stderr, "\n+++++++++++ ");
						printString(newId);
						fprintf(stderr, "\tgap: %d, counts: %d, N: %d, unsignificant (less than %d counts)\n", _gapIndex, countIUPAC[newId], N, minOcc+1);
					}
				}
			}
			if(depth < IUPAC_chars){
				/* recursivly go to more degenerated IUPAC strings */
				createIUPAC_list_phase1_rec(bestPhaseList, listSize, minOcc, N, newId, i+1, depth);
			}
		}
	}
}

void ElongCandidates::createIUPAC_list_phase2_rec(std::list<std::pair<uint32_t, double> >& bestPhaseList,
		int& listSize, int minOcc, int N, uint32_t id, MotifRegion& motifRegion, int pos, int depth){

	std::list<std::pair<uint32_t, double> >::iterator it;
	int endRange = Global::posSet->max_leng - (coding_characters+gaps[_gapIndex][coding_characters-2]) +1;
	region r;
	r.max=0;r.startRegion=0;r.endRegion=0;
	uint32_t newId, tmp = id;

	/* increment number of IUPAC characters in current string */
	depth++;
	/* shift to pos position */
	tmp = tmp >> pos*4;

	/* every position which is to degenerate */
	int changePositions = coding_characters;
	if(type == PALINDROME || type == TANDEM){
		changePositions = (coding_characters + 1) / 2;
	}

	for (int i = pos; i < changePositions; i++) {
		const uint8_t *IUPAC_base = IUPACS[tmp & 15];
		tmp = tmp >> 4;

		/* null out ith character */
		uint32_t erased_id = id & eraseBits[i];
		if(type == PALINDROME){
			erased_id = erased_id & eraseBits[coding_characters-i-1];
		}else if(type == TANDEM){
			if(i > 3) erased_id = erased_id & eraseBits[i-3];
			else	  erased_id = erased_id & eraseBits[i+3];
		}

		/* every allowed degeneration */
		for(int j=1; j<4; j++){
			/* add new iupac character at ith position */
			int base = IUPAC_base[j];
			newId = erased_id | (base << 4*i);
			if(type == PALINDROME){
				newId = newId | (IUPAC_palindrome[base] << 4*(coding_characters-i-1));
			}else if(type == TANDEM){
				if(i>3) newId = newId | (base << 4*(i-3));
				else    newId = newId | (base << 4*(i+3));
			}


			startPosList& current = collectStartPosHash[newId];

			if(current.getPhase2()){
				/* get region of motif enrichment */
				if(Global::usePositionalProbs){
					r = motifRegion.getRegion(current.getStartPosCount(), current.size());
				}else{
					r.set = 0;
					r.max = 1;
					r.startRegion = 1;
					r.endRegion = endRange;
				}
				//printString(newId);
				double set_p = getLog_Pvalue(newId, current, r.startRegion, r.endRegion, Global::multipleOccurrence);
				//cerr << ": " << set_p << "\t";

				if(listSize < (int)NUMBER_PHASE2_DEGENERATIONS){
					for(it = bestPhaseList.begin(); it != bestPhaseList.end(); it++){
						if(it->second > set_p)break;
					}
					bestPhaseList.insert(it, std::pair<uint32_t, double>(newId, set_p));
					listSize++;
					/* if motif has better pval than worst motif in best motif list, drop worst element and include new one
					 * at correct position */
				}else if((int)NUMBER_PHASE2_DEGENERATIONS == 1 || set_p < bestPhaseList.back().second){
					if(bestPhaseList.size() > 0) bestPhaseList.pop_back();
					for(it = bestPhaseList.begin(); it != bestPhaseList.end(); it++){
						if(it->second > set_p)break;
					}
					bestPhaseList.insert(it, std::pair<uint32_t, double>(newId, set_p));
				}
			}
			if(depth < IUPAC_characters){
				/* recursivly go to more degenerated IUPAC strings */
				createIUPAC_list_phase2_rec(bestPhaseList, listSize, minOcc, N, newId, motifRegion, i+1, depth);
			}
		}
	}
}


int ElongCandidates::setLog_Pvalues(){
	/* only calculate pValue if motif occurs in at least 3 sequences or 1% of all sequences */
	int minOcc = std::max(3, nbSequences / 100);

	const uint8_t* gap;
	if(type == FIVEMERS) gap = gaps[_gapIndex];
	else gap = gaps_palindrome[_gapIndex];
	int N = totalPositions - nbSequences * (coding_characters + gap[coding_characters-2] - 1);
	if(type == PALINDROME && Global::revcomp)
		N = totalPositions/2 - nbSequences * (coding_characters + gap[coding_characters-2] - 1);

	int maxSeqLength = Global::posSet->max_leng;

	int count = 0;

	unsigned int ids;
	if(type == FIVEMERS) ids = POW_4_coding;
	else if(type == PALINDROME || type == TANDEM) ids = POW_4_3;
	else {cerr << "motif type error !" << endl; exit(-1);}

	for(uint32_t i=0; i< ids; i++){
		uint32_t tmp = i;
		uint32_t newId = 0;
		if(type == FIVEMERS){
			for(int j=0; j<coding_characters; j++){
				int base = tmp&3; tmp = tmp >> 2;
				newId = (newId << 4) + base;
			}
		}else if(type == PALINDROME){
			for(int j=0; j<coding_characters/2; j++){
				int base = tmp&3; tmp = tmp >> 2;
				newId += (base << (coding_characters-j-1)*4);
				newId += ((IUPAC_palindrome[base]) << j*4);
			}
		}else if(type == TANDEM){
			for(int j=0; j<coding_characters/2; j++){
				int base = tmp&3; tmp = tmp >> 2;
				newId += (base << 4*j);
				newId += (base << 4*(j+3));
			}
		}

		std::list<std::pair<uint32_t, double> > bestPhaseList;

		double p = kmerProbs[newId];
		/* calculate pValue for finding countIUPAC[id] number of instances in set */
		double set_p = calculateOrderStatisticsPvalue(countIUPAC[newId], N, p) + LOG_Bonferonni;

		for(std::list<std::pair<uint8_t,uint32_t> >::iterator it = TEST_ID.begin(); it != TEST_ID.end(); it++){
			if(_gapIndex == (*it).first && newId == (*it).second){
				fprintf(stderr, "\n+++++++++++ ");
				printString(newId);
				fprintf(stderr, "\tgap: %d, counts: %d, N: %d, p: %f\t, E-Value: %e\n\n", _gapIndex, countIUPAC[newId], N, kmerProbs[newId], exp(set_p));
			}
		}

		bestPhaseList.push_back(std::pair<uint32_t, double>(newId, set_p));

		int listSize = 1;
		//cerr << endl; printString(newId); cerr << "(" << set_p << ") : ";
		createIUPAC_list_phase1_rec(bestPhaseList, listSize, minOcc, N, newId, 0, 0);

		/* add 10 best iupac degenerations to next round */
		for(std::list<std::pair<uint32_t, double> >::iterator it = bestPhaseList.begin(); it != bestPhaseList.end(); it++){
			uint32_t id = it->first;

			if(!collectStartPosHash[id].getPhase2()){
				collectStartPosHash[id].initialize(countIUPAC[id], N, maxSeqLength, Global::usePositionalProbs);
				collectStartPosHash[id].setPhase2();
				count++;
			}
		}
	}

	for(unsigned int id=0; id < POW_4_coding; id++){
		int tmp = id;
		uint32_t newId = 0;
		for(int j=0; j<coding_characters; j++){
			int base = tmp&3; tmp = tmp >> 2;
			newId = (newId << 4) + base;
		}
		/* always store start positions for non degenerative kmers */
		if(!collectStartPosHash[newId].getPhase2() && countIUPAC[newId] > 0){
			collectStartPosHash[newId].initialize(countIUPAC[newId], N, maxSeqLength, Global::usePositionalProbs);
			collectStartPosHash[newId].setPhase2();
		}
	}

	return count;
}

void ElongCandidates::createStartPositionList(){
	const uint8_t *gap = gaps[_gapIndex];
	if(type == PALINDROME || type == TANDEM){
		gap = gaps_palindrome[_gapIndex];
	}

	/* every sequence */
	for (int i = 1; i <= nbSequences; i++) {
		unsigned char *sequence = Global::posSet->entity[i]->S[0];
		/* every startPos */
		//cerr << i << "/" << Global::posSet->entity[i]->n << "/" << Global::posSet->entity[i]->n-coding_characters-gap[coding_characters-2]+1 << endl;
		int seqLen = Global::posSet->entity[i]->n-coding_characters-gap[coding_characters-2]+1;
		if(type == PALINDROME && Global::revcomp) seqLen /= 2;
		for (int j = 1; j <= seqLen; j++) {
			int base = sequence[j]-1;
			if (base < 0) continue;
			uint32_t id = base;

			bool N_char = false;

			/* read sequence backwards to generate id as later transformation of id from base 4 to base 16 reverses the sequence */
			for (int l = 1; l<coding_characters; l++) {
				int base = sequence[j+l+gap[l-1]]-1;
				if (base < 0){ N_char = true; break; }
				id = (id << 4) + base;
			}

			/* if sequence had no N character, store start position */
			if(!N_char){
				collectStartPosHash[id].push_back(i, j);

			}
		}
	}
	/* store start position in all IUPAC sequences */
	setStartPos();
}

void ElongCandidates::setStartPos(){
	for(uint32_t id=0; id < POW_4_coding; id++){
		uint32_t tmp = id;
		uint32_t newId = 0;
		for(int j=0; j<coding_characters; j++){
			int base = tmp&3; tmp = tmp >> 2;
			newId = (newId << 4) + base;
		}
		/* if not a single instance exist, nothing to merge */
		if(!collectStartPosHash[newId].capacity()) continue;

		/* merge start Positions */
		mergeIUPACS_rec(collectStartPosHash, newId, collectStartPosHash[newId], 0, 0);

		/* if nondegenerated kmer not significant, remove it */
		if(!collectStartPosHash[newId].getPhase2()) collectStartPosHash[newId].reset();
	}
}

int ElongCandidates::filterElongationCandidates(elongList &elongatedList){
	/* create motif Region object */

	int minOcc = std::max(3, nbSequences / 100);
	int N = nbSequences;

	int startRange = 1;
	int endRange = Global::posSet->max_leng - (coding_characters+gaps[_gapIndex][coding_characters-2]) +1;
	if(type == PALINDROME && Global::revcomp){
		endRange = Global::posSet->max_leng / 2 - (coding_characters+gaps[_gapIndex][coding_characters-2]) +1;
	}

	MotifRegion motifRegion(startRange, endRange);

	/* initialize a bitwise kmer representation with 4 bits for a nucleotide and 2 bits for gaps */
	const uint8_t *gap = gaps[_gapIndex];
	if(type == PALINDROME || type == TANDEM) gap = gaps_palindrome[_gapIndex];


	Extension::graphVizOptions *vizopts = NULL;
	std::ofstream* grout = NULL;

	if(TEST_ID.size() > 0){
		std::string trackingFile = Global::outputDirectory;
		trackingFile += "/extGraph.dot";
		grout = new std::ofstream(trackingFile.c_str());

		/* extend motif */
		vizopts = new Extension::graphVizOptions;
		vizopts->showEdgesToLosers = true;
		vizopts->showEdgesToLookups = true;
		vizopts->stream = grout;
		*(vizopts->stream) << "digraph G {" << endl;
	}

	region r;
	r.set=0;r.max=0;r.startRegion=0;r.endRegion=0;
	bool* filteredIUPACs = (bool*)calloc(POW_16_coding, sizeof(bool));

	int count = 0;

	unsigned int ids;
	if(type == FIVEMERS) ids = POW_4_coding;
	else if(type == PALINDROME || type == TANDEM) ids = POW_4_3;
	else {cerr << "motif type error !" << endl; exit(-1);}

	for(uint32_t i=0; i< ids; i++){
		std::list<std::pair<uint32_t, double> > bestPhaseList;

		uint32_t tmp = i;
		uint32_t newId = 0;
		if(type == FIVEMERS){
			for(int j=0; j<coding_characters; j++){
				int base = tmp&3; tmp = tmp >> 2;
				newId = (newId << 4) + base;
			}
		}else if(type == PALINDROME){
			for(int j=0; j<coding_characters/2; j++){
				int base = tmp&3; tmp = tmp >> 2;
				newId += (base << (coding_characters-j-1)*4);
				newId += ((IUPAC_palindrome[base]) << j*4);
			}
		}else if(type == TANDEM){
			for(int j=0; j<coding_characters/2; j++){
				int base = tmp&3; tmp = tmp >> 2;
				newId += (base << 4*j);
				newId += (base << 4*(j+3));
			}
		}
		startPosList& current = collectStartPosHash[newId];

		int listSize = 0;
		double PVal=0;
		if(current.getPhase2()){
			/* get region of motif enrichment */
			if(Global::usePositionalProbs){
				r = motifRegion.getRegion(current.getStartPosCount(), current.size());
			}else{
				r.set = 0;
				r.max = 1;
				r.startRegion = 1;
				r.endRegion = endRange;
			}
			//cerr << endl; printString(newId);
			PVal = getLog_Pvalue(newId, current, r.startRegion, r.endRegion, Global::multipleOccurrence);
    		//cerr << ": " << PVal << endl;

			if(current.size() > 2) listSize = 1;
		}

		createIUPAC_list_phase2_rec(bestPhaseList, listSize, minOcc, N, newId, motifRegion, 0, 0);

		if(current.getPhase2() && current.size() > 2)
			bestPhaseList.push_back(std::pair<uint32_t, double>(newId, PVal));

		elongList bestPvalList;

		bool trackMotif = false;
		for(std::list<std::pair<uint32_t, double> >::iterator it = bestPhaseList.begin(); it != bestPhaseList.end(); it++){
			uint32_t id = it->first;
			if(!filteredIUPACs[id]){
				uint32_t tmpId = id;
				startPosList& filteredIUPAC = collectStartPosHash[id];

				/* get region of motif enrichment */
				if(Global::usePositionalProbs){
					r = motifRegion.getRegion(filteredIUPAC.getStartPosCount(), filteredIUPAC.size());
				}else{
					r.set = 0;
					r.max = 1;
					r.startRegion = 1;
					r.endRegion = endRange;
				}

				uint8_t b3 = (tmpId & 15); tmpId = (tmpId >> 4);
				int g2 = gap[coding_characters-2] - gap[coding_characters-3];
				uint8_t b2 = (tmpId & 15); tmpId = (tmpId >> 4);
				int g1 = gap[coding_characters-3] - gap[coding_characters-4];
				uint8_t b1 = (tmpId & 15); tmpId = (tmpId >> 4);

				/* generate first kmer as trimer with two gap */
				SmallKmer *kmer = new SmallKmer(g1, g2, b1, b2, b3);

				/* mutate the kmer to have to correct length */
				for(int j=coding_characters-4; j>0; j--){
					int g = gap[j] - gap[j-1];
					uint8_t b = (tmpId & 15);tmpId = (tmpId >> 4);
					kmer->mutate(-g-1, b);
				}
				int g = gap[0];
				uint8_t b = (tmpId & 15);tmpId = (tmpId >> 4);
				kmer->mutate(-g-1, b);

				/* create object with all information about the kmer */
				std::shared_ptr<Kmer> ptr(new Kmer(kmer));
				/* copy start Positions in needed data structure*/
				filteredIUPAC.copyStartPosToMatchContainer(ptr->seeds, type, kmer->length());

				/* set region */
				ptr->enrichment = r;

				Extension::graphVizOptions *vizopt = NULL;

				bool debugMotif = false;
				for(std::list<std::pair<uint8_t,uint32_t> >::iterator it = TEST_ID.begin(); it != TEST_ID.end(); it++){
					if(_gapIndex == (*it).first && id == (*it).second){
						vizopts->id_suffix = "___" + ptr->getKmer()->bestNucString();
						vizopt = vizopts;

						trackMotif = true;
						debugMotif = true;
					}
				}

				/* calculate pValue */
				Extension::calibrateMotif(*ptr, debugMotif);

				if(debugMotif){
					fprintf(stderr, "\n\n+++++++++++ before elongation: ");
					cerr << ptr->getKmer()->bestNucString();
					fprintf(stderr, "\tcount: %d, E-Value: %e, pVal Set: %e\n", ptr->setSize, ptr->p_pos, exp(ptr->p_set));
					for(MatchContainer::iterator it = ptr->seeds.begin(); it != ptr->seeds.end(); it++){
						for(int i=0; i<ptr->getKmer()->length(); i++){
							fprintf(stderr, "%c", AlphaChar(Global::posSet->entity[it->seq]->S[0][it->pos+i], Global::A));
						}
						fprintf(stderr, "(%d/%d) ", it->seq, it->pos);
					}
					fprintf(stderr, "\n");
				}

				/* extend motifs */
				elongList extendedList = Extension::getBestExtensions(elongList(1, ptr), extendedMotifsHash, vizopt, false);

				for(std::list<std::pair<uint8_t,uint32_t> >::iterator it = TEST_ID.begin(); it != TEST_ID.end(); it++){
					if(_gapIndex == (*it).first && id == (*it).second){
						fprintf(stderr, "+++++++++++ after elongation: ");
						if(extendedList.size() == 0) continue;

						std::shared_ptr<Kmer> ptr2 = extendedList.front();
						cerr << ptr2->getKmer()->bestNucString();
						fprintf(stderr, "\tE-Value: %e, pVal Set: %e\n\n", ptr2->p_pos, exp(ptr2->p_set));
						MatchContainer& seeds = ptr2->seeds;
						fprintf(stderr, "set size: %d\t", ptr2->setSize);
						for(MatchContainer::iterator it = seeds.begin(); it != seeds.end(); it++){
							fprintf(stderr, "%d/%d ", it->seq, it->pos);
						}
						fprintf(stderr, "\n");
						Global::trackedMotifs.insert(ptr2->getKmer()->bestNucString());
						Global::trackedElongation = true;
						//exit(-1);
					}
				}

				if(extendedList.size() > 0){
					double newPval = extendedList.front()->p_set;

					if(bestPvalList.size() < NUMBER_PHASE3_DEGENERATIONS || newPval < bestPvalList.back()->p_set){
						/* add list of extended motifs at the end of elongCandidates,
						 * if pValue < 1e10 and if it is in top MIN_IUPAC_EXTENSIONS */
						unsigned int counter = 0;
						for(elongList::iterator it = bestPvalList.begin(); it != bestPvalList.end(); it++, counter++){
							if(counter > NUMBER_PHASE3_DEGENERATIONS){
								while(it != bestPvalList.end()){
									it = bestPvalList.erase(it);
								}
							}else if(newPval < (*it)->p_set){
								bestPvalList.insert(it, extendedList.begin(), extendedList.end());
								break;
							}
						}if(counter < NUMBER_PHASE3_DEGENERATIONS){
							bestPvalList.insert(bestPvalList.end(), extendedList.begin(), extendedList.end());
						}
					}
				}
				filteredIUPACs[id] = true;
				count++;
			}
		}
		bool elongMotif = false;
		elongList::iterator it = bestPvalList.begin();
		while(it != bestPvalList.end()){
			if(trackMotif){
				if (Global::isTracked((*it)->getKmer()->bestNucString())) {
					elongMotif = true;
					if((*it)->seeds.begin() == (*it)->seeds.end()){
						fprintf(stderr, "tracked Motif within best IUPAC extensions, but has no seeds => filtered out\n");
					}
				}
			}
			(*it)->seeds.clear(); 		/* seeds do not have to be copied, they are received directly from the sequences later */
			elongatedList.push_back(*it);

			it = bestPvalList.erase(it);
		}
		if(trackMotif){
			if(!elongMotif)	fprintf(stderr, "tracked Motif not within best IUPAC extensions\n");
		}
	}
	free(filteredIUPACs);

	if(TEST_ID.size() > 0){
		*(vizopts->stream) << "}" << endl;
		grout->close();
		delete grout;
		delete vizopts;
	}
	return count;
}

double ElongCandidates::getLog_Pvalue(uint32_t id, startPosList& spl, int startRange, int endRange, bool mops){

	int N = nbSequences;
	int sequencesWithMotif = spl.countStartPosInRange(N, startRange, endRange, mops, motifsPerSequenceCount);
	if(sequencesWithMotif < 3) return 1e100;

	/* pVal for finding the motif at a random position */
	double p = kmerProbs[id];

	if(mops){
		const uint8_t* gap;
		if(type == FIVEMERS) gap = gaps[_gapIndex];
		else gap = gaps_palindrome[_gapIndex];

		N = totalPositions - nbSequences*(coding_characters+gap[coding_characters-2]-1);
		if(type == PALINDROME && Global::revcomp) N = totalPositions/2 - nbSequences * (coding_characters + gap[coding_characters-2] - 1);

		if(Global::usePositionalProbs) N = (endRange - startRange + 1) * nbSequences;

	}else{
		/* pVal for finding the motif in a random sequence (use average length of sequences with motif) */
		if(Global::usePositionalProbs) p = 1-pow(1-p, endRange-startRange+1);
		else if(type == PALINDROME && Global::revcomp){ p = 1-pow(1-p, avgLength/2); }
		else p = 1-pow(1-p, avgLength);

		N = nbSequences;
	}

	return calculateOrderStatisticsPvalue(sequencesWithMotif, N, p) + LOG_Bonferonni;
}

/* converts a gapless string into its id */
uint32_t ElongCandidates::StringToId(const char* s){
	uint32_t id = 0;
	int base = 0;
	for(unsigned int i=0; i<strlen(s); i++){
		switch(s[i]){
			case 'A': base = 0;break;case 'C': base = 1;break;case 'G': base = 2;break;	case 'T': base = 3;break;
			case 'M': base = 4;break;case 'R': base = 5;break;case 'W': base = 6;break;case 'S': base = 7;break;
			case 'Y': base = 8;break;case 'K': base = 9;break;case 'N': base = 10;break;
		}
		id = (id << 4) + base;
	}
	return id;
}

void ElongCandidates::setTrackedMotif(const std::unordered_set<std::string>& motifs){
	TEST_ID.clear();

	for (auto it = motifs.begin(); it != motifs.end(); ++it) {
		const std::string& motif = *it;
		if(motif.length() == 0){
			continue;
		}

		int trackedGaps[100];
		for(int i=0; i<100; i++)trackedGaps[i] = 0;

		std::string trackedMotif;
		int gapPos=0;

		trackedMotif.push_back(motif.at(0));
		for(unsigned int i=1; i<motif.length(); i++){
			char c = motif.at(i);
			if(c == 'N' || c == '.' || c == 'X' ){
				trackedGaps[gapPos]++;
				continue;
			}
			trackedMotif.push_back(c);
			gapPos++;
		}
		if(trackedMotif.length() != coding_characters){
			if(type == FIVEMERS){
				fprintf(stderr, "\n\ntrackedMotif has to have %d non gap positions for nonpalindromic case. motif %s has %d\n\n", coding_characters, motif.c_str(), (int)trackedMotif.length());
			}else{
				fprintf(stderr, "\n\ntrackedMotif has to have %d non gap positions for palindromic case. motif %s has %d\n\n", coding_characters, motif.c_str(), (int)trackedMotif.length());
			}
			return;
			//exit(-1);
		}
		for(int i=1; i<coding_characters-1; i++){
			trackedGaps[i] += trackedGaps[i-1];
		}

		int gapId = -1;
		for(int i=0; i<Global::GAPS; i++){
			const uint8_t *gap = gaps[i];
			if(type == PALINDROME || type == TANDEM) gap = gaps_palindrome[i];
			bool found = true;
			for(int j=0; j<coding_characters-1; j++){
				if(trackedGaps[j] != gap[j] ){
					found = false;
				}
			}
			if(found){
				gapId = i;
				break;
			}
		}

		//fprintf(stderr, "motif: %s\t", trackedMotif.c_str());
		//for(int i=0; i<coding_characters-1; i++) fprintf(stderr, "%d\t", trackedGaps[i]);

		if(gapId == -1){
			if(type == FIVEMERS){
				fprintf(stderr, "\n\ntrackedMotif %s has not sampled gap positions for non-palindromic case. Try a higher gapNumber with --gaps option\n\n", motif.c_str());
			}else{
				fprintf(stderr, "\n\ntrackedMotif %s has not sampled gap positions for palindromic case\n\n", motif.c_str());
			}
		}

		if(gapId != -1){
			fprintf(stderr, "\n\n"\
					"*********************************\n"\
					"* motif %s is tracked over time\n"\
					"* tracking of extension phase written to file %s/extGraph.dot\n"\
					"* translate it to ps with: dot -Tps %s/extGraph.dot -o extGraph.ps\n"\
					"**********************************\n\n",
					motif.c_str(), Global::outputDirectory, Global::outputDirectory);

			TEST_ID.push_back(std::pair<uint8_t,uint32_t>(static_cast<uint8_t>(gapId),StringToId(trackedMotif.c_str())));
		}
	}
}

void ElongCandidates::printString(uint32_t id) {
	const uint8_t *gap = gaps[_gapIndex];
	if(type == PALINDROME || type == TANDEM) gap = gaps_palindrome[_gapIndex];
	char IUPAC[max_pattern_width+1];
	int pos = coding_characters + gap[coding_characters - 2];
	IUPAC[pos--] = '\0';
	int base = (id & 15);
	IUPAC[pos--] = IUPAC_CHARS[base];
	for (int i = 1; i < coding_characters; i++) {
		int gapSize = gap[coding_characters-i-1];
		if(coding_characters-i-2 >= 0) gapSize -= gap[coding_characters-i-2];
		for (int m = 0; m < gapSize; m++) {
//			IUPAC[pos--] = Global::aa ? 'X' : 'N';
			IUPAC[pos--] = 'N';
		}
		id = (id >> 4);
		base = (id & 15);
		IUPAC[pos--] = IUPAC_CHARS[base];
	}
	fprintf(stderr, "%s\t", IUPAC);
}
