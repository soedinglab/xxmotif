/*
 * sequence.h
 *
 *  Created on: Oct 12, 2009
 *      Author: eckhart
 */

#ifndef MPRO_SEQUENCE_H_
#define MPRO_SEQUENCE_H_

#include <stdint.h>
#include <vector>

class Sequence {

public:
	Sequence(const std::vector<uint8_t> &v) :
		seq(v) {
	}

	size_t size() const {
		return seq.size();
	}

	uint8_t& operator[](const int i) {
		return seq[i];
	}

	uint8_t operator[](const int i) const {
		return seq[i];
	}

private:
	std::vector<uint8_t> seq;
};

#endif /* SEQUENCE_H_ */
