#include "branch_and_bound-inl.h"

std::ostream& operator<<(std::ostream &os, const pwm_kmer_t &k) {
	os << "[" << k.kmer << ":" << std::scientific << k.probability << "]";
	return os;
}
