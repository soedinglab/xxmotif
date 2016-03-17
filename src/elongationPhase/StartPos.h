#ifndef STARTPOS_H_
#define STARTPOS_H_

#include "../memoryPool/pool_alloc.h"
#include <list>
#include <stdint.h>

class StartPos {
public:
	StartPos(){};
	StartPos(const int32_t s, const int32_t p) :
		seq(s), pos(p){
	}

	int32_t seq;
	int32_t pos;

	bool operator== (const StartPos& other) const {
		return seq == other.seq && pos == other.pos;
	}

	static bool cmpStartPos(const StartPos& s1, const StartPos& s2) {
		if(s1.seq < s2.seq) return true;
		else if(s1.seq > s2.seq) return false;
		else{
			if(s1.pos <= s2.pos) return true;
			else return false;
		}
	}

private:

};

typedef std::list<StartPos, Pool_alloc<StartPos> > StartPosContainer;
//typedef std::list<StartPos> StartPosContainer;

#endif /* STARTPOS_H_ */
