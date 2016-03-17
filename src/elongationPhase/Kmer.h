#ifndef STATE_KMER_RESULT_H_
#define STATE_KMER_RESULT_H_

#include "Match.h"
#include "../SmallKmer.h"
#include "../UniversalKmer.h"
#include "../AbstractKmer.h"
#include "../nucleotides/motifRegion.h"
#include "../aminoacids/MProGlobal.h"

#include <iostream>
#include <memory>

class Kmer {

	/** Container for all relevant results concerning a single kmer */
private:
	AbstractKmer *kmer;
	void operator=(const Kmer&);
	Kmer();
public:

	Kmer(AbstractKmer *k) :
		kmer(k), p_pos(1), p_set(1), setSize(std::numeric_limits<int>::min()) {
	}

	~Kmer() {
		delete kmer;
	}
	/**
	 * Copy constructor
	 */
	Kmer(const Kmer &other) {
		kmer = other.getKmer()->clone();
		p_pos = other.p_pos;
		p_set = other.p_set;
		setSize = other.setSize;
		enrichment = other.enrichment;
		debugString = other.debugString;
		seeds = other.seeds;
	}

	Kmer(const AbstractKmer& otherKmer){
		if (otherKmer.numMatches() == SmallKmer::maxNumMatches) {
			kmer = new UniversalKmer(otherKmer);
		}else {
			kmer = otherKmer.clone();
		}
		p_pos = 1;
		p_set = 1;
		setSize = std::numeric_limits<int>::min();
	}

	AbstractKmer* getKmer() {
		return kmer;
	}

	const AbstractKmer* getKmer() const {
		return kmer;
	}

	void setKmer(AbstractKmer* k);

	/** p value of best match */
	double p_pos;
	/** p value of whole set */
	double p_set;
	/** optimal subset size as determined by order statistics */
	int setSize;
	/** all matches in positive set considered significant (p<=p_pos) */
	MatchContainer seeds;
	/* region of motif enrichment */
	region enrichment;
	/* for debug information to be shown upon output */
	std::string debugString;

	bool isDuplicate(const Kmer &other) const;
	bool operator==(const Kmer &other) const {
		return isDuplicate(other);
	}
	bool operator!=(const Kmer &other) const {
		return !(*this==other);
	}

	/** comparison by comparison of corrected p-values */
	bool operator<(const Kmer &other) const {
		return p_set < other.p_set;
	}
	struct lessPtrClass {
		bool operator()(const std::shared_ptr<Kmer>& lhs,
				const std::shared_ptr<Kmer>& rhs) const {
			return *lhs < *rhs;
		}
	};
	static bool lessPtr(const std::shared_ptr<Kmer> &s1,
			const std::shared_ptr<Kmer> &s2) {
		return s1->p_set < s2->p_set;
	}
	struct eqPtrClass {
		bool operator()(const std::shared_ptr<Kmer>& lhs,
				const std::shared_ptr<Kmer>& rhs) const {
			return lhs == rhs;
		}
	};
	static bool eqPtr(const std::shared_ptr<Kmer> &s1,
			const std::shared_ptr<Kmer> &s2) {
		return *s1 == *s2;
	}
};

inline void Kmer::setKmer(AbstractKmer* k) {
	delete kmer;
	kmer = k;
}

inline bool Kmer::isDuplicate(const Kmer &other) const {
	if (p_set != other.p_set) {
		return false;
	} else {
		MatchContainer::const_iterator m1 = seeds.begin();
		MatchContainer::const_iterator m2 = other.seeds.begin();
		while (m1 != seeds.end()) {
			if (m2 == other.seeds.end()) {
				return false;
			} else {
				if (*m1 != *m2) {
					return false;
				}
				++m1;
				++m2;
			}
		}
		if (m2 != other.seeds.end()) {
			return false;
		} else {
			return true;
		}
	}
}

std::ostream& operator<<(std::ostream &os, const Kmer &res);

#endif /* STATE_KMER_RESULT_H_ */
