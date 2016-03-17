#ifndef EXTENSIONTABLE_H_
#define EXTENSIONTABLE_H_

#include "../aminoacids/MProGlobal.h"

#include <cmath>
#include <iostream>
#include <list>
#include <valarray>

class ExtensionTable {

//	friend std::ostream& operator<<(std::ostream&, const ExtensionTable&);

public:

	typedef std::list<uint8_t> stateList_type;

	static void initialize(MProGlobal::SeqInfo_t &_S, const double p_thresh, const double seqThresh, const double seqThresh2 = 0);

	static void freeTable() {
		free(tablePointer);
		free(table);
		free(counter);
	}

	static void reset() {
		memset(table, 0, setSize * stateSize * sizeof(bool));
		memset(counter, 0, stateSize * sizeof(int));
	}

	static void count(const int seq, const int aa, bool mops);

	static void getExtensionStates(int set_size, stateList_type& extstates, bool Debug);


private:
	static bool initialized;
	static bool** tablePointer;
	static bool* table;
	static int* counter;
	static MProGlobal::SeqInfo_t *S;
	static std::vector<std::list<int> > sigStates;
	static int setSize;
	static int stateSize;
	static double pThreshold;
	static double seqThreshold;
	static double seqThreshold2;
};



inline void ExtensionTable::count(const int seq, const int aa, bool mops) {
	if(mops){
		for (std::list<int>::const_iterator it = sigStates[aa].begin(); it!=sigStates[aa].end(); ++it) {
			counter[*it] ++;
		}
	}else{
		bool* pointer = tablePointer[seq];
		for (std::list<int>::const_iterator it = sigStates[aa].begin(); it!=sigStates[aa].end(); ++it) {
			counter[*it] += (pointer[*it] ^ true);
			pointer[*it] = true;
		}
	}
}

//std::ostream& operator<<(std::ostream&, const ExtensionTable&);

#endif
