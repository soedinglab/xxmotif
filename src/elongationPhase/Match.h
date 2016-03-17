#ifndef MATCH_H_
#define MATCH_H_

#include "../memoryPool/pool_alloc.h"
#include <list>
#include <iomanip>
#include <iosfwd>
#include <ios>
#include <iostream>
#include <stdint.h>


/** representation of match in a set of sequences */
class Match {
public:
	Match(const int32_t s, const int32_t p, const float sc=1.0) :
		seq(s), pos(p), score(sc) {
	}

	int32_t seq;
	int32_t pos;
	float score;
	//std::string debugString;

	bool operator<(const Match &other) const {
		return score < other.score;
	}

	bool operator== (const Match& other) const {
		return seq == other.seq && pos == other.pos;
	}

	bool operator!=(const Match &other) const {
		return !(*this==other);
	}

	static bool cmpMatches(const Match& s1, const Match& s2) {
		if(s1.seq < s2.seq) return true;
		else if(s1.seq > s2.seq) return false;
		else{
			if(s1.pos <= s2.pos) return true;
			else return false;
		}
	}


private:
	Match(){};
};

inline std::ostream& operator<<(std::ostream &os, const Match &m) {
	std::iostream::fmtflags flags = os.flags();
	std::streamsize prec = os.precision();
	os << "[seq: " << std::setw(3) << m.seq << ", pos: " << std::setw(5)
			<< m.pos << ", score: " << std::scientific << std::setprecision(2)
			<< m.score << "]";// << " " << m.debugString;
	os.flags(flags);
	os << std::setprecision((int)prec);
	return os;
}

/** Container for all relevant results concerning a single kmer */
typedef std::list<Match, Pool_alloc<Match> > MatchContainer;

#endif /* MATCH_H_ */
