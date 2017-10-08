//
//  Molecule.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Molecule__
#define __vagabond__Molecule__

#include <stdio.h>
#include <vector>
#include "mat3x3.h"
#include "shared_ptrs.h"
#include "Sampler.h"
#include <string>

struct vec3;

class Molecule
{
public:
	void addModel(ModelPtr model);
	void addAtom(AtomPtr atom);

	void addToMap(FFTPtr fft, mat3x3 _real2frac);

	virtual void summary();
	virtual void tieAtomsUp() {};
	virtual void refine(CrystalPtr target, RefinementType rType);
	void tiedUpScattering(double *tied, double *all);
	virtual void makePDB(std::string filename);
	virtual void graph(std::string graphName) {};
	virtual void differenceGraphs(std::string graphName, CrystalPtr diffCryst) {};
	void resetInitialPositions();
	void setAnchors();
	
	long int atomCount()
	{
		return atoms.size();
	}

	AtomPtr atom(int i)
	{
		return atoms[i];
	}

	void setChainID(std::string chain)
	{
		_chainID = chain;
	}

	std::string getChainID()
	{
		return _chainID;
	}

	virtual std::string getClassName()
	{
		return "Molecule";
	}

private:
	std::vector<ModelPtr> models;
	std::vector<AtomPtr> atoms;


	std::string _chainID;
};

#endif /* defined(__vagabond__Molecule__) */
