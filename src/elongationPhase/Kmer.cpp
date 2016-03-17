#include "Kmer.h"
#include "elongationCandidates.h"
#include "../NullModel.h"
#include <cstddef>

std::ostream& operator<<(std::ostream &os, const Kmer &res) {
		//os << "operator<< undefined for Kmer in DNA mode" << endl;
		std::string ks = res.getKmer()->toString(0, ElongCandidates::IUPAC_CHARS);
		os << endl << "IUPAC string: " << ks << endl;
		os << "\tSeeds: ";
		size_t count = 0;
		for (MatchContainer::const_iterator m_it = res.seeds.begin(); m_it
				!= res.seeds.end(); ++m_it) {
			os << (count > 0 ? "\t       " : "") << *m_it << " (";
			for (size_t i = 0; i < ks.length(); ++i) {
				if (ks[i] != 'N') {
					os << AlphaChar(Global::posSet->entity[m_it->seq]->S[0][m_it->pos+i], Global::A);
				} else {
					os << ".";
				}
			}
			os << ")";
			if (count < res.seeds.size() - 1) {
				os << endl;
			}
			++count;
		}
		os << endl;
		os << "nbSeeds: " << res.seeds.size() << endl;
	return os;
}
