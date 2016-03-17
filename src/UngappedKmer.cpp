#include "UngappedKmer.h"

int UngappedKmer::bitsPerChar;

std::ostream& operator<<(std::ostream& os, const UngappedKmer &k) {
	os << k.toString();
	return os;
}
