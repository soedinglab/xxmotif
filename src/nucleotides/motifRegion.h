#ifndef MOTIF_REGION_H_
#define MOTIF_REGION_H_

#include "../Globals.h"
#include "../LogTable.h"
#include "../utils.h"

struct region{
	uint8_t set;
	int max;
	int startRegion;
	int endRegion;
};

class MotifRegion{
public:
	MotifRegion( int startRange, int endRange,
				 double enrichmentThreshold = 0.1,
				 double significanceThreshold = -LogTable::LOG_i[1000]);

	region getRegion(int* startPosList, int counts);

private:
	int findMaximum(int* startPosCounts);
	double calcRegionPval(int startRegion, int endRegion, int length, int x, int N);

	double LOG_Bonferonni;

	static const int _windowSize = 5;

	int _positions[_windowSize+1];
	const int _startRange;
	const int _endRange;
	const double _enrichmentThreshold;
	const double _significanceThreshold;
};

inline double MotifRegion::calcRegionPval(int startRegion, int endRegion, int length, int x, int N){
	double p = (endRegion-startRegion+1.0)/length;
	return calculateOrderStatisticsPvalue(x, N, p);
}

#endif /* MOTIF_REGION_H_ */
