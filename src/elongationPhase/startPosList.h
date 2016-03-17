#ifndef START_POS_LIST_H_
#define START_POS_LIST_H_

#include "StartPos.h"
#include "Match.h"
#include "../Globals.h"
#include <stdlib.h>
#include <string.h>
#include <iostream>

using std::ostream;

class startPosList{
public:
	startPosList(){
		_phase2_summary = NULL;
		_capacity=0;
		_count = 0;
		_phase2=false;
	}

	~startPosList()	{
		freeMemory();
	}

	struct phase2_summary{
		bool _enrichedRegion;
		int _maxSeqLength;
		int*  _startPosCounts;

		StartPos *_startPositions;
	};
	typedef phase2_summary *p_phase2_summary;

	startPosList& operator+=(startPosList &other);
	friend ostream& operator<< (ostream &os, const startPosList &spl);

	void initialize(int size, int max_number_of_seq, int max_seq_length, bool enriched_region);
	int* getStartPosCount() { return _phase2_summary->_startPosCounts;}

	//int size()const { return _phase2_summary->_startPositions.size();}
	uint32_t size()const { return _count;}
	uint32_t capacity()const { return _capacity;}
	void reset();
	void freeMemory();

	void setPhase2(){ _phase2 = true; }
	bool getPhase2()const { return _phase2; }

	void copyStartPosToMatchContainer(MatchContainer& seeds, motif_type type, int length);

	//int getSequence(int nb)const {return _phase2_summary->_startPositions[nb].seq; }
	//int getPos(int nb)const { return _phase2_summary->_startPositions[nb].pos; }

	//MatchContainer& getStartPositions() { return _phase2_summary->_startPositions; }

	void push_back(int seq, int pos);
	int countStartPosInRange(const int maxSeq, const int startRange, const int endRange, bool mops, int* motifsPerSequenceCount);

private:
	p_phase2_summary _phase2_summary;
	uint32_t _count;
	uint32_t _capacity;
	bool _phase2;
};

/* copy startpositions at the end of the current list
 * only if current list passed the first filter step
 * update one occurrence per sequence counter and startPosCounts */
inline startPosList& startPosList::operator+=(startPosList &other){
	if(getPhase2()){
		int size = other.size();
		StartPos* otherSP = other._phase2_summary->_startPositions;

		StartPos* startPositions = _phase2_summary->_startPositions;
		//MatchContainer& startPositions = _phase2_summary->_startPositions;

		int* startPosCounts = _phase2_summary->_startPosCounts;
		bool enrichedRegion = _phase2_summary->_enrichedRegion;

		//MatchContainer& otherSP = other._phase2_summary->_startPositions;
		//for(MatchContainer::iterator it = otherSP.begin(); it != otherSP.end(); it++){
		for(int i=0; i < size; i++){
			int seq = otherSP[i].seq;
			int pos = otherSP[i].pos;

			assert(_count+i < _capacity);
			startPositions[_count+i] = StartPos(static_cast<int32_t>(seq), static_cast<int32_t>(pos));
			//int seq = it->seq;
			//int pos = it->pos;
			//startPositions.push_back(Match(seq, pos));
			if(enrichedRegion) startPosCounts[pos]++;
		}
		_count += size;
	}
    return *this;
}

inline void startPosList::initialize(int size, int max_number_of_seq, int max_seq_length, bool enriched_region){
	if(_phase2_summary == NULL){
		_phase2_summary = new phase2_summary;
		_phase2_summary->_maxSeqLength = max_seq_length+1;
		_phase2_summary->_enrichedRegion = enriched_region;
		if(enriched_region){
			_phase2_summary->_startPosCounts = (int*) calloc(_phase2_summary->_maxSeqLength, sizeof(int));
		}
	}
	_phase2_summary->_startPositions = new StartPos[size];
	//fprintf(stderr, "allocate %d start positions\n", size);

	_capacity = size;
}

inline void startPosList::reset(){
	if(_phase2_summary != NULL){
		free(_phase2_summary->_startPositions);
		_phase2_summary->_startPositions = NULL;
		if(_phase2_summary->_enrichedRegion){
			memset(_phase2_summary->_startPosCounts, '\0', _phase2_summary->_maxSeqLength*sizeof(int));
		}
	}
	_count = 0;
	_capacity = 0;
	_phase2 = false;
}

inline void startPosList::push_back(int seq, int pos){
	//_phase2_summary->_startPositions.push_back(Match(seq, pos));
	//fprintf(stderr, "pb %d\n", _count);
	assert(_count < _capacity);

	_phase2_summary->_startPositions[_count++] = StartPos(static_cast<int32_t>(seq), static_cast<int32_t>(pos));
	if(_phase2_summary->_enrichedRegion)_phase2_summary->_startPosCounts[pos]++;
}

inline int startPosList::countStartPosInRange(const int maxSeq, const int startRange, const int endRange, bool mops, int* motifsPerSequenceCount){
	if(mops){
		memset(motifsPerSequenceCount, 0, (Global::posSet->nent+1)*sizeof(int));
		int count = 0;
		StartPos* sp = _phase2_summary->_startPositions;

		for(uint32_t i=0; i<_count; i++){
			if(!Global::usePositionalProbs || (sp[i].pos >= startRange && sp[i].pos <= endRange)){
				motifsPerSequenceCount[sp[i].seq]++;
				if(motifsPerSequenceCount[sp[i].seq] <= Global::maxMotifsPerSequence){
					count++;
				}
			}
		}
		return count;
	}else{
		int count = 0;
		bool* seqOcc = (bool*)calloc(maxSeq+1, sizeof(bool));
		StartPos* sp = _phase2_summary->_startPositions;

		for(uint32_t i=0; i<_count; i++){
			if(!seqOcc[sp[i].seq] && sp[i].pos >= startRange && sp[i].pos <= endRange){
				seqOcc[sp[i].seq] = true;
				count++;
			}
		}
		free(seqOcc);
		return count;
	}
}

inline void startPosList::copyStartPosToMatchContainer(MatchContainer& seeds, motif_type type, int length){
	StartPos* sp = _phase2_summary->_startPositions;
	for(uint32_t i=0; i<_count; i++){
		//fprintf(stderr, "normal: %d/%d\n", sp[i].seq, sp[i].pos);
		seeds.push_back(Match(sp[i].seq, sp[i].pos));
		if(Global::revcomp && type == PALINDROME){
			//fprintf(stderr, "palin: %d/%d\n", sp[i].seq, Global::posSet->entity[sp[i].seq]->n - sp[i].pos - length + 2);
			seeds.push_back(Match(sp[i].seq, static_cast<int32_t>(Global::posSet->entity[sp[i].seq]->n - sp[i].pos - length + 2)));
		}
	}
	seeds.sort(Match::cmpMatches);
}

inline void startPosList::freeMemory(){
	if(_phase2_summary != NULL){
		/* _phase2_summary->_startPositions is freed outside of class
		 * pointer to data structure is still needed */
		if(_phase2_summary->_startPositions != NULL){
			delete[] _phase2_summary->_startPositions;
		}
		if(_phase2_summary->_enrichedRegion){
			free(_phase2_summary->_startPosCounts);
		}
		delete _phase2_summary;
	}
	_phase2_summary=NULL;
}

inline ostream& operator<< (ostream &os, const startPosList &spl){
//	MatchContainer& startPositions = spl._phase2_summary->_startPositions;
//	for(MatchContainer::iterator it = startPositions.begin(); it != startPositions.end(); it++){
//		os << it->seq << "/" << it->pos << "\t";
//	}
	for(uint32_t i=0; i<spl._count; i++){
		int seq = spl._phase2_summary->_startPositions[i].seq;
		int pos = spl._phase2_summary->_startPositions[i].pos;

		os << seq << "/" << pos << "\t";
	}
	os << "\n";
	return os;
}


#endif /* START_POS_LIST_H_ */
