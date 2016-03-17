#ifndef _MPRO_GLOBAL_H
#define _MPRO_GLOBAL_H

#include "alphabet.h"
#include "sequence.h"
#include "StateLib.h"

class MProGlobal {
public:
	typedef struct {
		/** underlying alphabet of amino acids or whatever */
		Alphabet alphabet;
		/** profile states */
		StateLib states;
		/** negative set */
		std::vector<Sequence> _negSet;
		/** positive set */
		std::vector<Sequence> _posSet;
		/** max distance to search for extension column (from left/right end of kmer) */
		int D_MAX;
	} SeqInfo_t;

private:
	static MProGlobal *instance;
  MProGlobal();
  ~MProGlobal();
	SeqInfo_t S;

	/* helpers for initialization */

public:
	static std::string getFullResourceName(const std::string &name);
	static const SeqInfo_t& getS() {
	  if (instance==NULL) instance = new MProGlobal();
		return instance->S;
	}
};

#endif
