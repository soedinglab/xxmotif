#ifndef _BRANCH_AND_BOUND_H
#define _BRANCH_AND_BOUND_H

#include "NullModel.h"
#include "UngappedKmer.h"
#include "refinementPhase/Motif.h"
#include "nucleotides/DnaHash.h"
#include "memoryPool/pool_alloc.h"

#include <unordered_map>
#include <list>
#include <map>
#include <vector>
#include <string>

using std::list;
using std::vector;
using std::map;



typedef std::unordered_map<uint64_t, double> aa_hash_t;
typedef DnaHash dna_hash_t;

struct pwm_kmer_t {
	pwm_kmer_t() : probability(log(1)) {};
	pwm_kmer_t(const UngappedKmer &k, const float p=1.0f) : kmer(k), probability(p) {};
	UngappedKmer kmer;
	float probability;
};
typedef std::list<pwm_kmer_t, Pool_alloc<pwm_kmer_t> > prefix_list_t;

template <class HASH>
class Kmer_Generator_Reloaded {

private:
	Kmer_Generator_Reloaded();
	Kmer_Generator_Reloaded(const Kmer_Generator_Reloaded&); /* forbidden for singleton */
	Kmer_Generator_Reloaded& operator=(const Kmer_Generator_Reloaded&); /* forbidden for singleton */

public:
	typedef uint64_t id_type;
	typedef list<uint64_t, Pool_alloc<id_type> > id_list_type;
	typedef map<int, id_list_type, std::less<int>, Pool_alloc<std::pair<int, id_list_type > > > mapType;
	bool debug;

	static Kmer_Generator_Reloaded& getInstance() {
		static Kmer_Generator_Reloaded<HASH> instance;
		return instance;
	}
	~Kmer_Generator_Reloaded();
	void reinitialize(const AbstractKmer &kmer, const int splitAfter, double threshold_left, double threshold_right);
	void reinitialize(double const * const * const pwm, const motif_columns_type& motif_columns, int pos_set_size, int total_positions);
	double get_kmer_probability_matched_only(const std::string &kmer);
	double get_kmer_probability_matched_only(unsigned char* kmer);
	double get_kmer_probability_all_positions(unsigned char* kmer);
	const mapType& getScoreMap() const {
		return scoreMapArray;
	}
	int getIdShift() const { return ID_SHIFT; }
	int getChrMask() const { return CHR_MASK; }

private:

	void reinitialize_common();

	int ID_SHIFT;
	int CHR_MASK;
	int MAX_KMER_SIZE;
	int hashNumber;
	int totalLength;
	double** PWM;
	double** PWM_lin;
	int posSetSize;
	int totalPositions;
	std::vector<int> gaps;
	int gapOffset; /* 0 if in left or MAX_KMER_SIZE if in right part */
	double** swappedPWM;
	int **rowSwap;
	int *colSwap;
	double* S_max;
	int T;
	int counts;
	int part;
	bool fillProbabilityHash;
	int length;
	double realPval;
	double* bgLog;
	std::unordered_map<int, float> rightKmerCond; /** conditional probabilities for right part */
	std::unordered_map<int, float> rightKmerCondComp; /** conditional compositional probabilities for right part */
	struct score_vector{
		score_vector(const int s, const id_list_type& list) {
			score = s;
			id_list = list;
		}
		int score;
		id_list_type id_list;
		bool operator<(const score_vector &other) const{
			return score < other.score;
		}
	};
	typedef vector<score_vector> vectorType;
	vectorType scoreIdVector[2];
	mapType scoreMapArray;
	HASH pValueHash[2];
	HASH probabilityHash;

	//string leftConsKmer;
	prefix_list_t prefixList;

	void set_kmer_to_Pvalue_hash(int start, int stop);
	void create_similar_kmers_rec(int j, double S_j, unsigned char* idx);
	double calculate_kmer_probability(const id_type id);
	int get_Kmer_score(unsigned char* kmer);
	void printSwapPWM(int **pwm);
	void printSwappedPWM(double **pwm, const int start, const int stop);
	template <class Alloc> static unsigned int binarySearch(const vector<score_vector, Alloc>& a, int S);

	float getRightKmerConditional(const UngappedKmer& k);
	float getRightKmerConditionalCompositional(const UngappedKmer& k);
  	void recalculatePrefixProbabilities(prefix_list_t &kmers);

};

std::ostream& operator<<(std::ostream &os, const pwm_kmer_t &k);

#endif
