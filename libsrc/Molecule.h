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
#include "AtomGroup.h"
#include <string>

struct vec3;

class Molecule : public AtomGroup
{
public:
	void addModel(ModelPtr model);

	void addToMap(FFTPtr fft, mat3x3 _real2frac);

	virtual void summary();
	virtual void tieAtomsUp() {};
	virtual void refine(CrystalPtr target, RefinementType rType);
	void tiedUpScattering(double *tied, double *all);
	virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal);
	virtual void graph(std::string graphName) {};
	virtual void differenceGraphs(std::string graphName, CrystalPtr diffCryst) {};
	void resetInitialPositions();
	void setAnchors();

	virtual void reportParameters();

	void setChainID(std::string chain)
	{
		_chainID = chain;
	}

	std::string getChainID()
	{
		return _chainID;
	}

	std::vector<vec3> getCentroidOffsets()
	{
		return _centroidOffsets;
	}

	std::vector<mat3x3> getRotationCorrections()
	{
		return _rotations;
	}

	std::vector<vec3> getRotationCentres()
	{
		return _centroids;
	}
	virtual std::string getClassName()
	{
		return "Molecule";
	}

	MoleculePtr shared_from_this()
	{
		AtomGroupPtr groupPtr = AtomGroup::shared_from_this();
		return ToMoleculePtr(groupPtr);
	}
protected:
	std::vector<vec3> _centroidOffsets;
	std::vector<vec3> _centroids; // after offset correction
	std::vector<mat3x3> _rotations;

private:
	std::vector<ModelPtr> models;

	std::string _chainID;

};

#endif /* defined(__vagabond__Molecule__) */
