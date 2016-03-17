#include "AlignStream.h"

#ifndef __CUSTOM_MULT_FASTA_H__
#define __CUSTOM_MULT_FASTA_H__

class CustMultFasta: public AlignStream {
private:

public:
	CustMultFasta() :
		AlignStream() {
	}
	CustMultFasta(const AlignStream& a) :
		AlignStream(a) {
	}
	~CustMultFasta() {
	}
	void read(std::istream &s);
	std::ostream & print(std::ostream &s) const;
	CustMultFasta& operator=(const AlignStream& rhs) {
		this->assign(rhs.begin(), rhs.end());
		return *this;
	}
};

#endif
