#ifndef DNAHASH_H_
#define DNAHASH_H_

class DnaHash {

private:

	double* val;
	const int num_elements;

public:
	DnaHash() : num_elements(static_cast<int>(pow(4, std::max(8, (Global::maxMatchPositions + 1) / 2)))) {
		val = (double*)calloc(num_elements,sizeof(double));
	}

	~DnaHash() {
		free(val);
	}

	double operator[] (const uint64_t index) const {
		return val[index];
	}

	double& operator[](const uint64_t index) {
		return val[index];
	}

	size_t size() {
		return num_elements;
	}

	void clear() {
		memset(val, 0, num_elements*sizeof(double));
	}

};

#endif /* DNAHASH_H_ */
