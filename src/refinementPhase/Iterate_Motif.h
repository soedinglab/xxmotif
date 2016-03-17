#ifndef ITERATE_MOTIF_H_
#define ITERATE_MOTIF_H_

#include "Motif.h"
#include "Sorted_Sites.h"
#include "StartPosUpdater-inl.h"
#include "../pValCalculation.h"
#include "../branch_and_bound-inl.h"
#include "../nucleotides/motifRegion.h"

class IterateMotif{
public:
	IterateMotif();
	~IterateMotif();

	void optimizeMotif(Motif* motif, bool keepPosSetSize, bool plus);
	//bool checkPalindrome(Motif* motif, double pseudo);
	bool maximizeMotifLength(Motif* motif, bool lengthCorrection, bool palin);

private:
	void fillElongationList(list<int>& elongColumns, Motif* motif, bool lengthCorrection, int offset);
	bool fillTmpColumn(int mutation, motif_columns_type& tmpColumns, Motif* motif);
	bool fillTmpColumn_noGaps(int front, int back, motif_columns_type& tmpColumns, Motif* motif);
	void printPWM(double** pwm, motif_columns_type& motifColumn, ostream &os) const;

	double** elongPWM;
	double** bestElongPWM;
	double skipRoundThreshold;
	int minMatchPositions;
};

inline void IterateMotif::fillElongationList(list<int>& elongColumns, Motif* motif, bool lengthCorrection, int offset){
	int firstMotifColumn = motif->getFirstMotifColumn();
	int lastMotifColumn = motif->getLastMotifColumn();

	elongColumns.push_back(1000);

	/* motif length has not converged, try different lengths */
	if(lengthCorrection){
		for(int i = firstMotifColumn - offset; i <= lastMotifColumn + offset; i++){
			elongColumns.push_back(i);
		}
	}
}

inline bool IterateMotif::fillTmpColumn_noGaps(int front, int back, motif_columns_type& tmpColumns, Motif* motif){

	int firstMotifColumn = motif->getFirstMotifColumn();
	int lastMotifColumn = motif->getLastMotifColumn();

	tmpColumns = motif->getMotifColumns();

	/* use old columns */
	if(front == 0 && back == 0){
	// change positions add the beginning
	}else if(front != 0 && back == 0){
		//change positions at the beginning
		if(front < 0){
			for(int i=firstMotifColumn-1; i>firstMotifColumn-1+front; i--){
				tmpColumns.push_front(i);
			}
			// remove positions at the beginning
		}else if(front > 0){
			for(int i=0; i<front; i++){
				tmpColumns.pop_front();
			}
		}
	// change positions at the end
	}else if(front == 0 && back != 0){
		// remove positions at the end
		if(back < 0){
			for(int i=0; i>back; i--){
				tmpColumns.pop_back();
			}
		// add positions at the end
		}else if (back > 0){
			for(int i=lastMotifColumn+1; i<lastMotifColumn+1+back; i++){
				tmpColumns.push_back(i);
			}
		}
	// shift motif to the left
	}else if(front == -1 && back == -1){
		tmpColumns.push_front(firstMotifColumn-1);
		tmpColumns.pop_back();
	// shift motif to the right
	}else if(front == 1 && back == 1){
		tmpColumns.pop_front();
		tmpColumns.push_back(lastMotifColumn+1);
	// skip round
	}else{
		return true;
	}
	/* check that motif is not shorter than minimum length */
	if((int)tmpColumns.size() < Global::minMatchPositions) return true;
	/* check that motif is not longer than maximum length */
	if((int)tmpColumns.size() > Global::maxMatchPositions) return true;

	return false;
}


inline bool IterateMotif::fillTmpColumn(int mutation, motif_columns_type& tmpColumns, Motif* motif){
	int firstMotifColumn = motif->getFirstMotifColumn();

	tmpColumns = motif->getMotifColumns();

	/* use old columns */
	if(mutation == 1000){
	}else if(mutation > PWM_LENGTH || mutation < 1) {
		return true;
	/* current column is before first match column */
	}else if(mutation < firstMotifColumn){
		tmpColumns.push_front(mutation);
	}else{
		motif_columns_type::iterator it_column = tmpColumns.begin();
		for( ;it_column != tmpColumns.end(); it_column++){
			/* if current column is a match column => remove it */
			if(mutation == *it_column){
				for(int a = 1; a <= nAlpha(Global::A); a++){
					if(motif->getPWM()[mutation][a] > skipRoundThreshold){
						return true;
					}
				}

				it_column = tmpColumns.erase(it_column);
				break;
			}
			/* if current column is not in list of match columns => insert it into list */
			else if(mutation < *it_column){
				tmpColumns.insert(it_column, mutation);
				break;
			}
		}
		/* if current column is behind last match column */
		if(it_column == tmpColumns.end() && mutation != motif->getLastMotifColumn()){
			tmpColumns.push_back(mutation);
		}
	}
	/* check that motif is not shorter than minimum length */
	if((int)tmpColumns.size() < Global::minMatchPositions) return true;
	/* check that motif is not longer than maximum length */
	if((int)tmpColumns.size() > Global::maxMatchPositions) return true;

	return false;
}

#endif /* ITERATE_MOTIF_ */
