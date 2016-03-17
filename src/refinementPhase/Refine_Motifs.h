#ifndef MOTIF_REFINEMENT_H_
#define MOTIF_REFINEMENT_H_

#include "MotifContainer.h"

class MotifRefinement{
public:
	MotifRefinement(MotifContainer& startModels) : _models(startModels){}

	void start();

	MotifContainer& getMotifs(){ return _models; }

private:
	void trimMatrices(){};

	MotifContainer& _models;
};

#endif /* MOTIF_REFINEMENT_H_ */
