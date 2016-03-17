#ifndef EM_H
#define EM_H

#include "../refinementPhase/MotifContainer.h"

#include "hoUtils.h"
#include "../seqFormat/Alignment.h"

class EM{
public:

	EM( MotifContainer& startModels ) : _models( startModels ){}

	MotifContainer& getModels(){
		return _models;
	}

	void go();

	void score( ss_type sequences, bool logProbs=false );
	void evaluatePWMs( ss_type sequences, bool logProbs=false );

private:

	MotifContainer& _models;
};

#endif /* EM_H */
