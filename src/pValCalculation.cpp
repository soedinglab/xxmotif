#include "pValCalculation.h"

PVal_Calculator::PVal_Calculator(){
	positionalProb = (double*)calloc(Global::posSet->max_leng+1, sizeof(double));
	sum_i = (double*)calloc(Global::posSet->max_leng+1, sizeof(double));
	for(int i=0; i<=Global::posSet->max_leng; i++){
		sum_i[i] = 0;
		for(int j=i+1;j<=Global::posSet->max_leng;j++){
			sum_i[i] += 1.0/j;
		}
	}

	tmpArray = (int**)malloc((MAX_CONSERVATION_LENGTH+1)*sizeof(int*));
	tmpArray2 = (int**)malloc((MAX_CONSERVATION_LENGTH+1)*sizeof(int*));
	for(int i=0; i<=MAX_CONSERVATION_LENGTH; i++){
		tmpArray[i] = (int*)calloc(Global::posSet->max_MultSeq, sizeof(int));
		tmpArray2[i] = (int*)calloc(Global::posSet->max_MultSeq, sizeof(int));
	}
	pConsCorrection = (double*)malloc(PWM_LENGTH*sizeof(double));
	for(int i=0; i<PWM_LENGTH; i++)pConsCorrection[i] = pow(Global::consCorrection, i);

	pOverrepCorrection = (double*)malloc(PWM_LENGTH*sizeof(double));
	for(int i=0; i<PWM_LENGTH; i++)pOverrepCorrection[i] = pow(Global::overrepCorrection, i);

	result_ali_free = (uint8_t*) malloc((Global::posSet->max_leng+1) * sizeof(uint8_t));

}

PVal_Calculator::~PVal_Calculator(){
	free(positionalProb);
	free(sum_i);

	for(int i=0; i<=MAX_CONSERVATION_LENGTH; i++){
		free(tmpArray[i]);
		free(tmpArray2[i]);
	}
	free(tmpArray);
	free(tmpArray2);
	free(pConsCorrection);
	free(pOverrepCorrection);
	free(result_ali_free);
}

void printI(uint32_t index){
	int A = index % 16;	index /= 16;
	int C = index % 16;	index /= 16;
	int G = index % 16;	index /= 16;
	int T = index % 16;
	fprintf(stderr, "A: %d, C: %d, G: %d, T: %d\n", A, C, G, T);
}

/* compare function for sorting sortedSites_type arrays */
inline int compare_uint8_t (const void *a, const void *b)
{
  int diff = (int)(*(uint8_t*)b) - (int)(*(uint8_t*)a);
  if(diff < 0.0) return 1;
  else if(diff > 0.0) return -1;
  return 0;
}

double PVal_Calculator::calculatePvalCons(int seq, int pos, int firstMotifColumn, const motif_columns_type& sortedColumns){
	e_type sequence = Global::posSet->entity[seq];
	const unsigned short multSeqNb = sequence->mseq;
	int seqLength = sequence->n;
	float pVal_ali_based = 1;
	float pVal_ali_free = 1;

	if(multSeqNb == 1){	return pVal_ali_based; }

	/* get multSeqNb in the motif area */
	int realMultSeqNb = 1;
	for(; realMultSeqNb < multSeqNb; realMultSeqNb++){
		motif_columns_type::const_iterator it = sortedColumns.begin();
		for(; it != sortedColumns.end(); it++){
			if(sequence->S[realMultSeqNb][pos + *it - firstMotifColumn] == 0){
				break;
			}
		}
		if(it != sortedColumns.end()) break;
	}

	bool debug = false;
//	if(seq < 100) debug = true;

	if(debug){
		fprintf(stderr, "sortedColumns: ");
		for(motif_columns_type::const_iterator it = sortedColumns.begin(); it != sortedColumns.end(); it++){
			fprintf(stderr, "%d\t", *it);
		}
		fprintf(stderr, "\n");
	}

	/***
	 * Alignment based pVal
	 ***/
	if(realMultSeqNb > 1){
		uint32_t index = 0;
		int mutations = 0, maxMutations = 0, k = 0;
		for(motif_columns_type::const_iterator it = sortedColumns.begin(); it != sortedColumns.end(); it++){
			int l = pos + *it - firstMotifColumn;

			int base = sequence->S[0][l];
			index += mutIndex[base-1];
			for(uint8_t m=1; m<realMultSeqNb; m++){
				// don't count sequences where the part is not alignable
				if(sequence->S[m][l] != base) mutations++;
				//fprintf(stderr, "l: %d, m: %d, base: %d vs %d\n", *it-firstMotifColumn, m, base, sequence->S[m][l]);
			}

			maxMutations += realMultSeqNb-1;
			if(++k == MAX_CONSERVATION_LENGTH) break;
		}

		mutations = std::min(std::min((int)(maxMutations*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS), mutations);
		pVal_ali_based = Global::conservationProbs[index][realMultSeqNb-1][mutations];
		//pVal_ali_based = Global::conservationProbs[k][realMultSeqNb-1][mutations];
		if(debug)	fprintf(stderr, "\nali based mutations %d: pVal: %f, realMultSeqNb: %d, seq: %d, pos: %d\n", mutations, pVal_ali_based, realMultSeqNb, seq, pos);

	}

	/***
	 * Alignment free pVal
	 ***/
	if(Global::useAliFree){
		uint32_t index = 0;
		int k = 0, maxSize = 0;

		/* initialize query that all bases are correct */
		uint8_t q[80];
		memset(q, '\0', 80*sizeof(uint8_t));
		for(motif_columns_type::const_iterator it = sortedColumns.begin(); it != sortedColumns.end(); it++){
			int base = sequence->S[0][pos + *it - firstMotifColumn];
			/* works just with motifs shorter or equal 16
			 * TODO build striped version */
			if(*it-firstMotifColumn < 16){
				/* insert mismatches into query profile for Alignment free counting */
				for (int a=0; a<=4; a++){
					if(a == base) continue;
					q[a*16 + *it-firstMotifColumn] = 1;
				}
				index += mutIndex[base-1];
				if(*it-firstMotifColumn > maxSize) maxSize = *it-firstMotifColumn;
				if(++k == MAX_CONSERVATION_LENGTH) break;
			}
		}

		int mutations = 0;
		int region = 50;
		int ali_free_start = std::max(1, pos-region);
		int ali_free_end = std::min(seqLength, pos+region);
		int ali_free_length = ali_free_end - ali_free_start + 1;
		int ali_free_pos = pos-region > 0 ? region : pos;

		for(uint8_t m=1; m<multSeqNb; m++){
			//memset(result_ali_free, '\0', PWM_LENGTH*sizeof(uint8_t));
			int muts = sse_count_mismatches_gaps(q, &sequence->S[m][ali_free_start], ali_free_length, maxSize, ali_free_pos);
			if(Global::revcomp){
				int ali_free_start_rev = std::max(1, seqLength-pos-region);
				int ali_free_end_rev = std::min(seqLength, seqLength-pos+region);
				int ali_free_length_rev = ali_free_end_rev - ali_free_start_rev + 1;
				int ali_free_pos_rev = seqLength-pos-region > 0 ? region : seqLength-pos;
				int muts_rev = sse_count_mismatches_gaps(q, &sequence->S[m][ali_free_start_rev], ali_free_length_rev, maxSize, ali_free_pos_rev);
				//fprintf(stderr, "muts: %d, muts_rev: %d\n", muts, muts_rev);
				muts = std::min(muts, muts_rev);
			}
			mutations += muts;
		}
		mutations = std::min(std::min((int)(k*(multSeqNb-1)*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS), mutations);

		pVal_ali_free = Global::alignmentFreeProbs[index][multSeqNb-1][mutations];

		if(debug){
			fprintf(stderr, "AliFree: ");printI(index);
			fprintf(stderr, "seq: %d, pos: %d, multSeq: %d, mutations: %d => pVal: %f\n", seq, pos, multSeqNb-1, mutations, pVal_ali_free);
		}

		//if(pVal_ali_free > 0.5) pVal_ali_free = 1;


		if(debug){
			fprintf(stderr, "aliBased: => %d/%d => pVal: %f\n", seq, pos, pVal_ali_based);
			fprintf(stderr, "aliFree: => %d/%d => pVal: %f\n", seq, pos, pVal_ali_free);
		}
		//return pVal_ali_based*pVal_ali_free;
		//return std::min(pVal_ali_based, pVal_ali_free);
		return pVal_ali_free;
	}else{
		return pVal_ali_based;
	}

}


/* precalculate positional probabilities */
void PVal_Calculator::initPositionalProbs(int motifLength, int firstMotifColumn, region r, int offset){
	int length = Global::posSet->max_leng - motifLength + 1;

	if(Global::revcomp){
		length = (Global::posSet->max_leng / 2 ) - motifLength + 1;
	}

	int motifStart = r.startRegion + firstMotifColumn - 1 + offset;
	int motifEnd = r.endRegion + firstMotifColumn - 1 + offset;
	int startRegion = std::max(1, motifStart);
	int endRegion = std::min(motifEnd, length);

	//cout << "length: " << length << "\tstartRegion: " << startRegion << "\tendRegion: " << endRegion << "\n";

	for(int i=1;i<=length;i++){
		int dist=0;
		if(i < startRegion) dist = startRegion - i;
		else if(i > endRegion) dist = i - endRegion;
		int min = std::max(1, startRegion - dist);
		int max = std::min(length, endRegion + dist);

		positionalProb[i] = (max-min+1.0)/length;
		//cerr << i << " " << positionalProb[i] << "\t";
	}


	if(Global::revcomp)	{
		// mirror the positional pvalue distribution on second half of sequence
		int start = Global::posSet->max_leng / 2 + 2;
		int oldLength = length;
		int length = Global::posSet->max_leng - motifLength + 1;

		int spaceToMiddleStart = Global::posSet->max_leng/2 - (motifStart + motifLength - 1);
		int spaceToMiddleEnd = Global::posSet->max_leng/2 - (motifEnd + motifLength - 1);

		startRegion = std::max(start, start + spaceToMiddleEnd );
		endRegion = std::min(length, start + spaceToMiddleStart );

		//cout << "\nlength: " << length << "\tstartRegion: " << startRegion << "\tendRegion: " << endRegion << "\n";

		for(int i=start;i<=length;i++){
			int dist=0;
			if(i < startRegion) dist = startRegion - i;
			else if(i > endRegion) dist = i - endRegion;
			int min = std::max(start, startRegion - dist);
			int max = std::min(length, endRegion + dist);


			positionalProb[i] = (max-min+1.0)/(oldLength);
			// cerr << i << " " << positionalProb[i] << "\t";
		}
	}
	//cerr << endl;
}

motif_columns_type PVal_Calculator::getSortedColumns(double** pwm, motif_columns_type& columns){
	motif_columns_type sorted_columns;
	double maxValue[PWM_LENGTH+1];
	double threshold = log(0.5);
	//double* bgLog = (Global::negSet==NULL) ? Global::posBg_log : Global::negBg_log;
	for(int i=1; i<=PWM_LENGTH; i++){
		//maxValue[i] = pwm[i][1] - bgLog[1];
		maxValue[i] = pwm[i][1];
		for(int j=2; j<=nAlpha(Global::A); j++){
			//if(pwm[i][j] - bgLog[j] > maxValue[i]) maxValue[i] = pwm[i][j] - bgLog[j];
			if(pwm[i][j] > maxValue[i]) maxValue[i] = pwm[i][j];
		}
	}

	for(unsigned int i = 0; i<columns.size(); i++){
		double max = -std::numeric_limits<double>::max();
		int bestColumn = 0;
		for(motif_columns_type::iterator it = columns.begin(); it != columns.end(); it++){
			if(maxValue[*it] > max){
				max = maxValue[*it];
				bestColumn = *it;
			}
		}
		//cerr << bestColumn << "(" << maxValue[bestColumn] << ")\t";
		//if(maxValue[bestColumn] > LogTable::LOG_2){
		if((sorted_columns.size() < 1 || maxValue[bestColumn] > threshold) && sorted_columns.size() < MAX_CONSERVATION_LENGTH){
			sorted_columns.push_back(bestColumn);
		}
		maxValue[bestColumn] = -std::numeric_limits<float>::max();
	}
	//cerr << endl;

	return sorted_columns;
}

PVal_Calculator& PVal_Calculator::getInstance()
{
  static PVal_Calculator instance;
  return instance;
}
