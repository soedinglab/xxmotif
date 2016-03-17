#ifndef __CONSERVATION_H__
#define __CONSERVATION_H__

#include "../Globals.h"
#include <stdint.h>

/********************************* PUBLIC *******************************/
void setConservationProbs();

/********************************* PRIVATE *******************************/
float*** allocateArray(ss_type set);
void printIndex(uint32_t index);

inline float*** allocateArray(ss_type set){
	float*** probs = (float***)calloc(((int)pow(2,16)+1), sizeof(float**));
	/* initialize Array that can index at most A + C + G + T = MAX_CONSERVATION_LENGTH kmers
	 * e.g. A=10, C=0, G=0, T=0 */

	for(int A=0; A<=MAX_CONSERVATION_LENGTH; A++){
		for(int C=0; C<=MAX_CONSERVATION_LENGTH; C++){
			for(int G=0; G<=MAX_CONSERVATION_LENGTH; G++){
				for(int T=0; T<=MAX_CONSERVATION_LENGTH; T++){
					int length = A+C+G+T;
					if(length>MAX_CONSERVATION_LENGTH) continue;
					uint32_t index = A*mutIndex[0] + C*mutIndex[1] + G*mutIndex[2] + T*mutIndex[3];
					probs[index] = (float**)malloc(set->max_MultSeq * sizeof(float*));
					length = std::min(length, MAX_CONSERVATION_LENGTH);
					for(int m=1; m<set->max_MultSeq; m++){
						probs[index][m] = (float*)calloc(std::min((int)(length*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS) + 1, sizeof(float));
					}
				}
			}
		}
	}
	return probs;
}

inline int sse_count_mismatches(const uint8_t* query_profile,
								const uint8_t* db_sequence,
								const int      dbseq_length,
								char* results, int pos)
{

return 0;
}

#endif
