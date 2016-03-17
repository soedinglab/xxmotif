#include "ThresholdChecker.h"

#include <cstdio>
#include <cstdlib>

ThresholdChecker::ThresholdChecker(const std::string &ts) : spec(ts) {
	size_t pos;
	std::string s(spec);
	if (0 == (pos = spec.find("<="))) {
		mode = LE;
		s.replace(0, 2, "");
	} else if (0 == (pos = spec.find("<"))) {
		mode = L;
		s.replace(0, 1, "");
	} else {
		fprintf(stderr, "Illegal threshold specification: %s\n", spec.c_str());
		throw(std::string("Illegal threshold specification"));
	}
	thresh = atof(s.c_str());
}

bool ThresholdChecker::satisfies(const double value) const {
	switch (mode) {
	case L:
		return value < thresh;
	case LE:
		return value <= thresh;
	default:
		fprintf(stderr, "ERROR: ThresholdChecker not initialized\n");
		throw(std::string("ThresholdChecker not initialized"));
	}
}
