#ifndef UNGAPPEDKMER_H_
#define UNGAPPEDKMER_H_

#include "AbstractKmer.h"
#include <limits>

class UngappedKmer: public AbstractKmer {

public:

	typedef std::list<UngappedKmer, Pool_alloc<UngappedKmer> >
			ungapped_kmer_list_t;

	UngappedKmer() {
		_length = 0;
		_numMatches = 0;
		id.val.num = 0;
		id.isNumeric = true;
	}

	UngappedKmer(const unsigned char * const s, const int len) {
		int res = *s;
		for (int i = 1; i < len; ++i) {
			res <<= bitsPerChar;
			res += *(s + i);
		}
		id.val.num = res;
		id.isNumeric = true;
		_length = len;
		_numMatches = len;
	}

	UngappedKmer(const std::string &s) {
		int res = AlphaCode((unsigned char)s[0], Global::A);
		int firstMatch = res == 0 ? std::numeric_limits<int>::max() : 0;
		for (int i = 1; i < static_cast<int>(s.length()); ++i) {
			res <<= bitsPerChar;
			res += AlphaCode((unsigned char)s[i], Global::A);
			if (firstMatch == std::numeric_limits<int>::max() && res != 0) {
				firstMatch = i;
			}
		}
		id.val.num = res;
		id.isNumeric = true;
		_length = std::max(1, (int) s.length() - firstMatch);
		_numMatches = _length;
	}



	operator int() const {
		return static_cast<int>(id.val.num);
	}

	inline uint8_t charAt(const int n) const {
		static const int bits = (1 << bitsPerChar) - 1;
		const int shift = (length() - 1 - n) * bitsPerChar;
		return static_cast<uint8_t>((id.val.num & (bits << shift)) >> shift);
	}

	inline void mutate(const int offset, const uint8_t res) {
		if (offset >= length()) {
			const int shift = (offset - length() + 1) * bitsPerChar;
			id.val.num <<= shift;
			id.val.num |= res;
			_length = offset + 1;
			_numMatches = offset + 1;
		} else if (offset >= 0) {
			static const int bits = (1 << bitsPerChar) - 1;
			const int shift = (length() - 1 - offset) * bitsPerChar;
			id.val.num &= ~(bits << shift);
			id.val.num |= res << shift;
		} else {
			const int shift = (length() - 1 - offset) * bitsPerChar;
			id.val.num |= res << shift;
			_length += -offset;
			_numMatches = _length;
		}
	}

	inline UngappedKmer mutatedCopy(const int offset, const uint8_t res) const {
		UngappedKmer copy(*this);
		copy.mutate(offset, res);
		return copy;
	}

	inline int gapsAfter(const int n) const {
		return 0;
	}

	motif_columns_type getMotifColumns() const {
		motif_columns_type motif_columns;
		for (int i = 0; i < numMatches(); i++) {
			motif_columns.push_back(i);
		}
		return motif_columns;
	}

	UngappedKmer* clone() const {
		return new UngappedKmer(*this);
	}

	inline bool operator<=(const AbstractKmer &other) const {
		return id.val.num <= other.getId().val.num;
	}

	UngappedKmer& operator++() {
		bool incremented = false;
		int digit = 0;
		while (!incremented) {
			const uint8_t res = charAt(length() - 1 - digit);
			if (res < nAlpha(Global::A)) {
				mutate(length() - 1 - digit, static_cast<uint8_t>(res + 1));
				incremented = true;
			} else {
				mutate(length() - 1 - digit, 0);
				++digit;
			}
		}
		return *this;
	}

	UngappedKmer& next_noWildcard() {
		bool incremented = false;
		int digit = 0;
		while (!incremented) {
			const uint8_t res = charAt(length() - 1 - digit);
			if (res < nAlpha(Global::A)) {
				mutate(length() - 1 - digit, static_cast<uint8_t>(res + 1));
				incremented = true;
			} else {
				mutate(length() - 1 - digit, 1);
				++digit;
			}
		}
		return *this;
	}

	/**
	 * Concatenation operator.
	 */
	const UngappedKmer operator*(const UngappedKmer &k2) const {
		uint64_t new_id = id.val.num;
		new_id <<= bitsPerChar * k2.length();
		new_id += k2.id.val.num;
		return UngappedKmer(new_id, _length+k2.length());
	}
	const UngappedKmer operator*(const char c) const {
		uint64_t new_id = id.val.num;
		new_id <<= bitsPerChar;
		new_id += AlphaCode(c, Global::A);
		return UngappedKmer(new_id, _length+1);
	}

	inline UngappedKmer subKmer(const int start, const int stop) {
		const int sub_len = stop - start + 1;
		assert(start>=0 && stop-length()<0);
		int bits = (1 << (sub_len * bitsPerChar)) - 1;
		bits <<= (length() - 1 - stop) * bitsPerChar;
		return UngappedKmer((id.val.num & bits) >> (length() - 1 - stop) * bitsPerChar, stop-start+1);
	}

	inline uint64_t getLastPosCombined(const int lastPos, const UngappedKmer& k){
		uint64_t bits = (1 << (lastPos * bitsPerChar)) -1;
		return ((id.val.num & bits) << k.length()*bitsPerChar) + k;
	}

	inline int lastWildcardPosition() const {
		static const int bits = (1 << bitsPerChar) - 1;
		for( int l=0; l <= length()-1; l++){
			if((id.val.num & (bits << l*bitsPerChar)) == 0){
				return length()-1-l;
			}
		}
		return undef;
	}

	ungapped_kmer_list_t expansions() const {
		ungapped_kmer_list_t exp(1, *this);
		expansions_rec(exp);
		return exp;
	}

	static inline UngappedKmer smallestWithLength(const int len) {
		assert(len>0);
		return UngappedKmer(1 << (bitsPerChar * (len - 1)), len);
	}

	static inline UngappedKmer greatestWithLength(const int len) {
		if (!(len > 0)) {
			fprintf(stderr, "Length %d\n", len);
			//			char *dummy; *dummy = 'c'; // provoke segfault
		}
		assert(len>0);
		int l = len;
		uint64_t idx(0);
		while (l > 1) {
			idx += nAlpha(Global::A);
			idx <<= bitsPerChar;
			--l;
		}
		idx += nAlpha(Global::A);
		return UngappedKmer(idx, len);
	}

	static inline UngappedKmer smallestWithLength_noWildcard(const int len) {
		assert(len>0);
		int l = len;
		uint64_t idx(0);
		while (l > 1) {
			idx += 1;
			idx <<= bitsPerChar;
			--l;
		}
		idx += 1;
		return UngappedKmer(idx, len);
	}

	static inline UngappedKmer greatestWithLength_noWildcard(const int len) {
		return greatestWithLength(len);
	}

	static const int undef = -2147483647;

	std::string toString(const int numBestCols = 0, const char* alphabet = NULL) const {
		std::stringstream s;
		for (int pos = 0; pos < length(); ++pos) {
			const int res = charAt(pos);
			if (res > nAlpha(Global::A)) {
				return std::string();
			} else {
				s << AlphaChar(res, Global::A);
			}
		}
		return s.str();
	}

	static void setNumberOfBits(const int n) {
		bitsPerChar = n;
		//fprintf(stderr, "Number of bits (UngappedKmer): %d\n", bitsPerChar);
	}

	static bool containsDollar(const UngappedKmer &kmer) {
		static const uint8_t dollar_index = AlphaCode('$', Global::A);
		for (int i = 0; i < kmer.length(); ++i) {
			if (kmer.charAt(i) == dollar_index) {
				return true;
			}
		}
		return false;
	}

private:

	static int bitsPerChar;

	UngappedKmer(const uint64_t index, const int l) {

		_length = std::max(1, l);
		_numMatches = _length;
		id.val.num = index;
		id.isNumeric = true;
	}

	inline static void expansions_rec(ungapped_kmer_list_t &exp) {
		ungapped_kmer_list_t::iterator it = exp.begin();
		while (it != exp.end()) {
			const int lx = (*it).lastWildcardPosition();
			if (lx == undef) {
				++it;
			} else {
				for (uint8_t res = 1; res <= nAlpha(Global::A); ++res) {
					UngappedKmer k(*it);
					k.mutate(lx, res);
					ungapped_kmer_list_t e2(1, k);
					expansions_rec(e2);
					exp.insert(it, e2.begin(), e2.end());
				}
				it = exp.erase(it);
			}
		}
	}

};

std::ostream& operator<<(std::ostream& os, const UngappedKmer &k);

#endif /* UNGAPPEDKMER_H_ */
