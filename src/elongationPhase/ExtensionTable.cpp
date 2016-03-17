#include "ExtensionTable.h"

#include "../utils.h"
#include <iomanip>
using std::endl;

bool ExtensionTable::initialized;
bool** ExtensionTable::tablePointer;
bool* ExtensionTable::table;
int* ExtensionTable::counter;
MProGlobal::SeqInfo_t *ExtensionTable::S;
std::vector<std::list<int> > ExtensionTable::sigStates;
int ExtensionTable::setSize;
int ExtensionTable::stateSize;
double ExtensionTable::pThreshold;
double ExtensionTable::seqThreshold;
double ExtensionTable::seqThreshold2;

void ExtensionTable::initialize(MProGlobal::SeqInfo_t &_S, const double p_thresh, const double seqThresh,
		const double seqThresh2) {
	assert(!initialized);
	reset();
	S = &_S;
	pThreshold = p_thresh;
	seqThreshold = seqThresh;
	seqThreshold2 = seqThresh2;

	sigStates.resize(S->alphabet.size());

	stateSize = static_cast<int>(S->states.size());

	for (int i = 0; i < static_cast<int>(S->alphabet.size()); ++i) {
		if (i == S->alphabet.wildcardIndex()) {
			continue;
		}
		sigStates[i].clear();
		for (int s = 0; s < stateSize; ++s) {
				if (S->states[s][i] >= pThreshold) {
					sigStates[i].push_back(s);
				}
		}
	}

	setSize = static_cast<int>(S->_posSet.size());
	tablePointer = (bool**)calloc(setSize, sizeof(bool*));
	table = (bool*)calloc(setSize*stateSize, sizeof(bool));
	counter = (int*)calloc(stateSize, sizeof(int));
	bool* ptr = table;
	for (int i = 0; i < setSize; ++i, ptr += stateSize) {
		tablePointer[i] = ptr;
	}
	initialized = true;
}

void ExtensionTable::getExtensionStates(int set_size, stateList_type& extStates, bool Debug) {
	double s_thresh = std::max((seqThreshold * set_size), Global::minCoverage);
	for (uint8_t s = 0; s < stateSize; ++s) {
		const int count = counter[s];

		if (s == 4 ) {
			s_thresh = (seqThreshold2 * set_size);
		}
		if (count >= s_thresh) {
			if(Debug)cerr << "state: " << (int)s << ": " << count << "\tsetSize: " << set_size << "\tthresh: " << s_thresh << endl;
			//if(count > set_size) exit(-1);
			extStates.push_back(s);
		}
		//if(Debug)cerr << "state: " << s << ": " << count << "\tsetSize: " << set_size << "\tthresh: " << s_thresh << endl;
	}
}
