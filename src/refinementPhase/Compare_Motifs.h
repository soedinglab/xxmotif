#ifndef MOTIF_COMPARISON_H_
#define MOTIF_COMPARISON_H_

#include "Motif.h"
#include "../Globals.h"

class MotifComparison{
public:
	static int calcBestMergingDist(Motif* model1, Motif* model2, int &bestRevcomp);
	static void mergePWMwithPWM(Motif* &model1, Motif* model2, int merge, int mergeRevcomp);
	static void removeOverlappingMotifOccurrencies(Motif* model1, Motif* model2, int merge);

	static bool DEBUG;
private:
	static int fillMergedModel(Motif* model1, Motif* model2, Motif* mergedModel, int offset, int tolerance);
	static int fillMergedModelRevcomp(Motif* model1, Motif* model2, Motif* mergedModel, int offset, int tolerance);
	static double calc_intermotif_dist(double probPWM1[31][21], double probPWM2[31][21], Motif* model1, Motif* model2, int revcomp, int offset, double max_PWM_diversity);
};


/* compare function for sorting sortedSites_type arrays */
inline int compare_doubles (const void *a, const void *b)
{
  double diff = *(double*)b - *(double*)a;
  if(diff > 0.0) return 1;
  else if(diff < 0.0) return -1;
  return 0;
}

inline int MotifComparison::fillMergedModel(Motif* model1, Motif* model2, Motif* mergedModel, int offset, int tolerance){
	StartPosContainer& merged_sites  = mergedModel->getStartPosList();
	StartPosContainer::iterator it_merged = merged_sites.begin();

	StartPosContainer& sites2  = model2->getStartPosList();
	StartPosContainer::const_iterator it2 = sites2.begin();

	int redundantMotifs = 0;

	while(it2 != sites2.end()){
		const int32_t newStart = static_cast<int32_t>(it2->pos + offset); //+model2->getFirstMotifColumn()-merge;
		const int32_t seq = it2->seq;

		//fprintf(stderr, "seq: %d, startPos: %d vs seq: %d, startPos: %d\n", it2->seq, newStart, it_merged->seq, it_merged->pos);
		if(it_merged == merged_sites.end() || seq < it_merged->seq || (seq == it_merged->seq && newStart < it_merged->pos-tolerance)){
			it_merged = merged_sites.insert(it_merged, StartPos(seq, newStart));
		}else{
			while( seq > it_merged->seq || (seq == it_merged->seq && newStart > it_merged->pos+tolerance) ){
				it_merged++;
				if(it_merged == merged_sites.end())break;
			}
			if(it_merged != merged_sites.end() && seq == it_merged->seq && abs(newStart - it_merged->pos) <= tolerance){
				//fprintf(stderr, "motif %d / %d is redundant\n", seq, newStart);
				redundantMotifs++;
			}else{
				it_merged = merged_sites.insert(it_merged, StartPos(seq, newStart));
			}
		}
		it_merged++;
		it2++;
	}
	return redundantMotifs;
}

inline int MotifComparison::fillMergedModelRevcomp(Motif* model1, Motif* model2, Motif* mergedModel, int offset, int tolerance){
	StartPosContainer& merged_sites  = mergedModel->getStartPosList();
	StartPosContainer::iterator it_merged = merged_sites.begin();

	StartPosContainer& sites2  = model2->getStartPosList();
	StartPosContainer::const_iterator it2 = sites2.begin();

	int redundantMotifs = 0;

	while(it2 != sites2.end()){
		const int32_t seq = it2->seq;
		StartPosContainer::const_iterator start = it2;

		// go to the lowest start position in sequence
		while(it2 != sites2.end() && it2->seq == seq)it2++;

		// fill mergedMotifs list
		do{
			it2--;
			/* go to the position which is the reverse complement of the first position of motif 1 */
			int32_t newStart = static_cast<int32_t>(it2->pos + model1->getFirstMotifColumn() - 1 + offset + model2->getMotifLength() - 1);
			/* calculate position on the reverse strand and go back to the start position of the PWM */
			newStart = static_cast<int32_t>(Global::posSet->entity[it2->seq]->n - newStart + 1 - model1->getFirstMotifColumn() + 1);

			if(DEBUG)fprintf(stderr, "seq: %d, startPos: %d vs seq: %d, startPos: %d\n",	it2->seq, newStart, it_merged->seq, it_merged->pos);
			if(it_merged == merged_sites.end() || seq < it_merged->seq || (seq == it_merged->seq && newStart < it_merged->pos-tolerance)){
				it_merged = merged_sites.insert(it_merged, StartPos(it2->seq, newStart));
			}else{
				while( seq > it_merged->seq || (seq == it_merged->seq && newStart > it_merged->pos+tolerance) ){
					it_merged++;
					if(it_merged == merged_sites.end()) break;
				}
				if(it_merged != merged_sites.end() && seq == it_merged->seq && abs(newStart-it_merged->pos) <= tolerance){
					//fprintf(stderr, "motif %d / %d is redundant\n", seq, newStart);
					redundantMotifs++;
				}else{
					it_merged = merged_sites.insert(it_merged, StartPos(seq, newStart));
				}
			}
			it_merged++;
		}while(it2 != start);

		// go to next sequence
		while(it2 != sites2.end() && it2->seq == seq) it2++;
	}
	return redundantMotifs;
}

inline double MotifComparison::calc_intermotif_dist(double probPWM1[31][21], double probPWM2[31][21], Motif* model1, Motif* model2,
		int revcomp, int offset, double max_PWM_diversity){
	int minOverlap = 6;

	double* bgLog = Global::negSet == NULL ? Global::posBg_log : Global::negBg_log;
	double bestDist = 1;
	int borderPos = 3;

	for(int pos1_offset = 0; pos1_offset < model1->getLastMotifColumn(); pos1_offset++){
		double dist = 0;
		int overlap = 0;
		double entropies[PWM_LENGTH];

		for(int pos1=std::max(1, model1->getFirstMotifColumn()-borderPos) + pos1_offset; pos1 <= std::min(PWM_LENGTH, model1->getLastMotifColumn()+borderPos); pos1++){
			int pos2 = pos1 + offset;

			// if reverse complement, look from the back
			if(revcomp == 1) {
				pos2 = model1->getFirstMotifColumn() + offset + model2->getMotifLength() - 1 - (pos1 - model1->getFirstMotifColumn());
			}

			if(pos2 < model2->getFirstMotifColumn()  || pos2 > model2->getLastMotifColumn()){ continue; }
			//if(DEBUG)fprintf(stderr, "\npos1: %d, pos2: %d\n", pos1, pos2);
			double sum = 0;
			double entropy = 0;
			for(int nuc1=1; nuc1 <= 4; nuc1++){
			  // if reverse complement, choose complementary nucleotide
			  int nuc2 = nuc1;
			  if(revcomp == 1){ nuc2 = 5-nuc2; }
			  double exp2 = probPWM2[pos2][nuc2];
			  double diff = probPWM1[pos1][nuc1] - exp2;

			  //if(DEBUG)fprintf(stderr, "%d: (%.2f, %.2f), diff: %.2f\t", nuc1, exp(model1->getPWM()[pos1][nuc1]), exp2, diff);

			  sum += diff*diff;


			  if(exp2 != 0){
				  entropy += exp2*(model2->getPWM()[pos2][nuc2]-bgLog[nuc2]) / LogTable::LOG_2;
			  }
			}
			entropies[overlap] = entropy;
			dist += sqrt(0.5 * sum);
			overlap++;

			if(overlap < minOverlap && (overlap != std::min(model1->getMotifLength(), model2->getMotifLength())) ){   continue;  }
			if(dist / overlap > max_PWM_diversity){	  continue;   }

				qsort(entropies, overlap, sizeof(double), compare_doubles);
				double avgEntropy = 0;
				int usedOverlap = std::min(overlap, minOverlap);
				for(int i=0; i< usedOverlap; i++){
					avgEntropy += entropies[i] / usedOverlap;
				}
				if(avgEntropy < 0.5) continue;

			if(dist/overlap < bestDist){
				bestDist = dist/overlap;
			}
		}
	}

	return bestDist;
}

#endif /* MOTIF_COMPARISON_H_ */
