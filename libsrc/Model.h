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
#include "Distributor.h"

// Anything which is capable of predicting electron positions.
//

class Model : public std::enable_shared_from_this<Model>, public Distributor
{
public:
	virtual FFTPtr getDistribution() = 0;

	virtual void addToMonomer(MonomerPtr monomer);
	virtual void addToMolecule(MoleculePtr molecule);

	virtual std::string getClassName() = 0;
	virtual vec3 getStaticPosition() = 0;

	void addDownstreamAtom(AtomPtr atom)
	{
		_downstreamAtoms.push_back(atom);
	}

	FFTPtr getZeroDistribution();

	std::vector<AtomWkr> _downstreamAtoms;
private:
	std::vector<AtomPtr> atoms;
};

#endif /* defined(__vagabond__Model__) */
