#include "UniversalKmer.h"

#include "AbstractKmer.h"
#include "Globals.h"
#include "SmallKmer.h"
#include "elongationPhase/Kmer.h"
#include "elongationPhase/Match.h"
#include "elongationPhase/elongationCandidates.h"
#include "memoryPool/pool_alloc.h"


#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdint.h>

const uint8_t UniversalKmer::charshift = 48;
const uint8_t UniversalKmer::gapChar = 252;

UniversalKmer::UniversalKmer(const int d, const uint8_t c1, const uint8_t c2) :
	cg() {
	cg.push_back(rep_type(c1, d));
	cg.push_back(rep_type(c2, 0));
	_length = 2 + d;
	_numMatches = 2;
	_strIndex = std::string(1, c1 + charshift);
	_strIndex += std::string(d, (char)gapChar);
	_strIndex += std::string(1, c2 + charshift);
	id.isNumeric = false;
	id.val.str = new char[_strIndex.length() + 1];
	strcpy(id.val.str, _strIndex.c_str());
}

UniversalKmer::UniversalKmer(const int g1, const int g2, const uint8_t c1,
		const uint8_t c2, const uint8_t c3) :
	cg() {
	cg.push_back(rep_type(c1, g1));
	cg.push_back(rep_type(c2, g2));
	cg.push_back(rep_type(c3, 0));
	_length = 3 + g1 + g2;
	_numMatches = 3;
	_strIndex = std::string(1, c1 + charshift);
	_strIndex += std::string(g1, (char)gapChar);
	_strIndex += std::string(1, c2 + charshift);
	_strIndex += std::string(g2, (char)gapChar);
	_strIndex += std::string(1, c3 + charshift);
	id.isNumeric = false;
	id.val.str = new char[_strIndex.length() + 1];
	strcpy(id.val.str, _strIndex.c_str());
}

UniversalKmer::UniversalKmer(const AbstractKmer& small) :
	cg(), _strIndex("") {
	_numMatches = small.numMatches();
	_length = small.length();
	for (int i = 0; i < small.numMatches(); ++i) {
		const rep_type elem(small.charAt(i),
				(i + 1 < small.numMatches()) ? small.gapsAfter(i) : 0);
		_strIndex += std::string(1, elem.chr + charshift);
		_strIndex += std::string(elem.gap, (char)gapChar);
		cg.push_back(elem);
	}
	id = AbstractKmer::id_type(_strIndex.c_str());
}

void UniversalKmer::mutate(const int offset, const uint8_t state) {
	if (offset < 0) { // append to front, fill with gaps
		cg.insert(cg.begin(), rep_type(state, -offset - 1));
		_length += -offset;
		_numMatches++;
		std::string pref(1, state + charshift);
		pref += std::string(-offset - 1, (char)gapChar);
		_strIndex.insert(0, pref);
	} else {
		int pos = 0;
		bool mutated = false;
		for (list_type::iterator it = cg.begin(); it != cg.end(); ++it) {
			if (offset == pos) { // mutate existing match position
				it->chr = state;
				mutated = true;
				_strIndex[pos] = static_cast<char>(state + charshift);
				break;
			} else if (pos + it->gap >= offset) { // insert here
				int oldgap = it->gap;
				it->gap = offset - pos - 1;
				_numMatches++;
				_strIndex.at(offset) = static_cast<char>(state + charshift);
				//_strIndex[offset] = state + charshift;
				list_type::iterator it2 = it;
				++it2;
				cg.insert(it2, rep_type(state, oldgap - it->gap - 1));
				mutated = true;
				break;
			}
			pos += 1 + it->gap;
		}
		if (!mutated) { // grow to the right
			cg.back().gap = offset - length();
			_strIndex += std::string(cg.back().gap, (char)gapChar);
			_strIndex += std::string(1, state + charshift);
			cg.push_back(rep_type(state, 0));
			_length = offset + 1;
			_numMatches++;
		}
	}
	delete[] id.val.str;
	id.val.str = new char[_strIndex.length() + 1];
	strcpy(id.val.str, _strIndex.c_str());
}

