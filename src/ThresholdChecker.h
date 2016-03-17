#ifndef THRESHOLDCHECKER_H_
#define THRESHOLDCHECKER_H_

#include <string>

class ThresholdChecker {

	enum Mode_t { L, LE };

private:
	Mode_t mode;
	double thresh;
	std::string spec;

public:
	ThresholdChecker() {};
	ThresholdChecker(const std::string &ts);
	bool satisfies(const double value) const;
	std::string getSpecification() const {
		return spec;
	}

};

#endif /* THRESHOLDCHECKER_H_ */
