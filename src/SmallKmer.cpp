#include "SmallKmer.h"

#include "alphabet.h"
#include "AbstractKmer.h"
#include "elongationPhase/Kmer.h"
#include "elongationPhase/Match.h"
#include "elongationPhase/elongationCandidates.h"
#include "memoryPool/pool_alloc.h"

#include "Globals.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <limits>
#include <math.h>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>

int SmallKmer::chrBits;
int SmallKmer::gapBits;
const uint64_t SmallKmer::allMask = std::numeric_limits<uint64_t>::max();
std::vector<uint64_t> SmallKmer::chrMask;
std::vector<uint64_t> SmallKmer::gapMask;
std::vector<uint64_t> SmallKmer::prefixMask;
std::vector<uint64_t> SmallKmer::suffixMask;

std::vector<int> SmallKmer::bitsBeforeChar;
std::vector<int> SmallKmer::bitsBeforeGap;
int SmallKmer::maxNumMatches;

bool SmallKmer::initialized;

std::string SmallKmer::convBase(uint64_t v, long base) {
	static const std::string digits = "0123456789abcdef";
	std::string result;
	if ((base < 2) || (base > 16)) {
		result = "Error: base out of range.";
	} else {
		do {
			result = digits[static_cast<int>(v % base)] + result;
			v /= base;
		} while (v);
	}
	return result;
}

void SmallKmer::init(int cbits, int gbits) {
	//assert(!initialized);
	chrBits = cbits;
	gapBits = gbits;
	maxNumMatches = 1 + (64 - chrBits) / (chrBits + gapBits);
	chrMask.resize(maxNumMatches);

	LOG(DEBUG3) << "chrBits: " << (int)chrBits << " gapBits: " <<  (int)gapBits  << " maxNumMatches: " << (int)maxNumMatches << endl;
	LOG(DEBUG3) << "             " << convBase(std::numeric_limits<uint64_t>::max(), 2) << endl;

	//	cout << "             " << convBase(std::numeric_limits<uint64_t>::max(), 2) << endl;

	const char fillchar = cout.fill();
	LOG(DEBUG3) << std::setfill('0');

	bitsBeforeChar.resize(maxNumMatches);
	for (uint8_t i = 0; i < chrMask.size(); ++i) {
		bitsBeforeChar[i] = i * (chrBits + gapBits);

		chrMask[i] = static_cast<uint64_t> (pow(2, chrBits) - 1) << bitsBeforeChar[i];
		LOG(DEBUG3) << "chrMask[" << std::setw(2) << (int) i << "]: " << std::setw(64) << convBase(chrMask[i], 2) << endl;
	}
	gapMask.resize(maxNumMatches - 1);
	bitsBeforeGap.resize(maxNumMatches - 1);
	for (int i = 0; i < static_cast<int>(gapMask.size()); ++i) {
		bitsBeforeGap[i] = ((i + 1) * chrBits + i * gapBits);
		gapMask[i] = static_cast<uint64_t> (pow(2, gapBits) - 1) << bitsBeforeGap[i];
		LOG(DEBUG3) << "gapMask[" << std::setw(2) << (int) i << "]: " << std::setw(64) << convBase(gapMask[i], 2) << endl;
	}
	prefixMask.resize(maxNumMatches);
	suffixMask.resize(maxNumMatches);
	prefixMask[0] = chrMask[0];
	for (uint8_t i = 1; i < prefixMask.size(); ++i) {
		prefixMask[i] = prefixMask[i - 1] | gapMask[i - 1] | chrMask[i];
	}
	for (uint8_t i = 0; (uint8_t) (i + 1) < suffixMask.size(); ++i) {
		suffixMask[i] = ~(prefixMask[i] | gapMask[i]);
	}
	suffixMask[suffixMask.size() - 1] = chrMask[chrMask.size() - 1];

	initialized = true;
	cout.fill(fillchar);
}

SmallKmer::SmallKmer(const int gap, const uint8_t c1, const uint8_t c2) {
	assert(initialized);
	_length = gap + 2;
	_numMatches = 2;
	id.isNumeric = true;
	id.val.num = allMask & ~chrMask[0];
	id.val.num += c2;
	id.val.num <<= gapBits;
	id.val.num += gap;
	id.val.num <<= chrBits;
	id.val.num += c1;
}


SmallKmer::SmallKmer(const int gap1, const int gap2, const uint8_t c1,const uint8_t c2, const uint8_t c3) {
	assert(initialized);
	_length = gap1 + gap2 + 3;
	_numMatches = 3;
	id.isNumeric = true;
	id.val.num = allMask & ~chrMask[0];
	id.val.num += c3;
	id.val.num <<= gapBits;
	id.val.num += gap2;
	id.val.num <<= chrBits;
	id.val.num += c2;
	id.val.num <<= gapBits;
	id.val.num += gap1;
	id.val.num <<= chrBits;
	id.val.num += c1;
}

void SmallKmer::mutate(const int offset, const uint8_t state) {
	if (offset < 0) { // grow left
		assert(_numMatches< maxNumMatches);
		assert(std::abs((float)offset)<=pow(2,gapBits));
		id.val.num <<= gapBits;
		id.val.num += -offset - 1;
		id.val.num <<= chrBits;
		id.val.num += state;
		_numMatches++;
		_length += -offset;
	} else if (offset >= length()) { // grow right
		assert(_numMatches<maxNumMatches);
		assert(offset-length()+1<=pow(2,gapBits));
		uint64_t n = state;
		n <<= gapBits;
		n += offset - length();
		n <<= bitsBeforeChar[numMatches()] - gapBits;
		n |= (allMask & ~(chrMask[numMatches()] | gapMask[numMatches() - 1]));
		id.val.num &= n;
		_numMatches++;
		_length += offset - length() + 1;
	} else {
		int pos = 0;
		for (uint8_t i = 0; i < numMatches(); ++i) {
			if (pos == offset) { // replace existing match position
				id.val.num &= ~chrMask[i];
				//						std::cerr << endl << stringId() << endl;
				uint64_t n = state;
				n <<= bitsBeforeChar[i]; //i * (chrBits + gapBits);
				id.val.num += n;
				break;
			} else {
				pos += 1 + gapsAfter(i);
				if (pos > offset) { // split gap
					assert(_numMatches<maxNumMatches);
					int oldgap = static_cast<int> ((id.val.num & gapMask[i])
							>> bitsBeforeGap[i]);
					int gap_after = pos - offset - 1;
					int gap_before = oldgap - gap_after - 1;
					uint64_t id_new = (id.val.num & suffixMask[i])
							>> bitsBeforeGap[i];
					id_new += gap_after;
					id_new <<= chrBits;
					id_new += state;
					id_new <<= gapBits;
					id_new += gap_before;
					id_new <<= bitsBeforeGap[i];
					id.val.num = id_new | (id.val.num & prefixMask[i]);
					_numMatches++;
					break;
				}
			}
		}
	}
}

std::string SmallKmer::asBinaryString() const {
	std::stringstream s;
	s << std::setfill('0');
	for (int i = numMatches() - 1; i >= 0; --i) {
		s << std::setw(chrBits) << convBase(charAt(i), 2);
		if (i > 0) {
			s << " " << std::setw(gapBits) << convBase(gapsAfter(i - 1), 2)
					<< " ";
		}
	}
	return s.str();
}
