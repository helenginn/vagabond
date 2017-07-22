//
//  Model.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Model__
#define __vagabond__Model__

#include "shared_ptrs.h"
#include <stdio.h>
#include <vector>


// Anything which is capable of predicting electron positions.
//

class Model : public std::enable_shared_from_this<Model>
{
public:
	virtual FFTPtr getDistribution(FFTPtr *reuse = NULL) = 0;

	virtual void addToMolecule(MoleculePtr molecule);

private:
	std::vector<AtomPtr> atoms;
};

#endif /* defined(__vagabond__Model__) */
