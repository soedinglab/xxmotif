#ifndef MOTIF_CONTAINER_H
#define MOTIF_CONTAINER_H

#include "Motif.h"
#include "Compare_Motifs.h"
#include "Iterate_Motif.h"
#include "../ThresholdChecker.h"

class MotifContainer{
	friend class Output;
public:
	MotifContainer(){};
	~MotifContainer();

	void add(Motif* m){ _startModels.push_back(m); }
	int getMotifNb() const { return static_cast<int>(_startModels.size()); }
	void clear();

	void update_Pval_and_sort(double pseudo);
	void updatePWM(double pseudo);
	void setBindingSites(const ThresholdChecker &pValThreshold);

	void merge(bool iterativePhase);
	void iterate(int maximizeMotifLength, bool lastIteration);

	void sortByPosPval();
	void sort();
	void sort_and_filter(double filterThreshold);

	void filter(const int _minimumMotifs, const int _maximumMotifs, const double _pValThreshold);
	void removeHomologousMotifs();
	void removeRedundantMotifs();

	/*
	 * filter by model number in ranking
	 */
	void filter( std::vector<int> nrModels );
	/*
	 * filter min. percentage of sequences containing a binding site instance
	 * in addition to
	 * min. number of models
	 * max. number of models
	 * max. p-value of models
	 * */
	void filter( const int minModels, const int maxModels, const double maxPvalue, const float minOccurrence );
	/*
	 * initialize higher-order models from XXmotif results
	 * */
	void initHoMotifsWithStartPos();
	/*
	 * calculate the percentage of positive sequences containing a binding site
	 * instance used as model-specific specificity factor in EM calculations
	 */
	void setSpecificityFactors();

//  /* Garbage motifs can seem very significant due to inaccuracies introduced by the split
//   * in the branch-and-bound calibration. These motifs should have drastically different
//   * P-values than their reverse complements, while true motifs should not.
//   * The filter finds motifs for which the difference between the regular motif
//   * and its reverse complement exceeds a threshold and removes them if desired.
//   *
//   * @params
//   * threshold: cutoff in log10 space
//   * remove: remove bad motifs or not
//   */
//  void revcompFilter();

	const list<Motif*>& getMotifs() const { return _startModels; }

private:

	list<Motif*> _startModels;
};

#endif /* MERGE_MOTIFS_H */
