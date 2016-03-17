#include "MProGlobal.h"

#include <iomanip>
#include <limits>
#include <string>
#include <vector>
#include <fstream>

MProGlobal* MProGlobal::instance;

MProGlobal::MProGlobal() {
	S.D_MAX = 7;
	S.alphabet = Alphabet(Global::A);
	S.states = StateLib(S.alphabet);
	/*
	 * convert ugly C structures to something more usable
	 */
	if (Global::negSet != NULL) {
		for (int i = 1; i <= Global::negSet->nent; i++) {
			const e_type &ent = Global::negSet->entity[i];
			std::vector<uint8_t> v((uint8_t*) (ent->S[0] + 1),
					(uint8_t*) (ent->S[0] + ent->n + 1));
			S._negSet.push_back(Sequence(v));
		}
	}
	for (int i = 1; i <= Global::posSet->nent; i++) {
		const e_type &ent = Global::posSet->entity[i];
		std::vector<uint8_t> v((uint8_t*) (ent->S[0] + 1),
				(uint8_t*) (ent->S[0] + ent->n + 1));
		S._posSet.push_back(Sequence(v));
	}
}

MProGlobal::~MProGlobal() {
  delete instance;
}

std::string MProGlobal::getFullResourceName(const std::string &name) {
	std::stringstream s;
	size_t pos = Global::argv0.find("Release/XXmotif");
	if ( pos == std::string::npos ) {
		pos = Global::argv0.find("Debug/XXmotif");
	}
	if ( pos == std::string::npos ) {
		pos = Global::argv0.find("GTest/XXmotif");
	}
	if ( pos != std::string::npos ) {
		std::string path(Global::argv0);
		path.erase(pos);
		s << path << "resources/mpro/" << name;
	} else {
		cerr << "ERROR: MPro resources not found, check your installation" << endl
				<< "argv0 is '" << Global::argv0 << "'" << endl;
		exit(1);
	}
	return s.str();
}
