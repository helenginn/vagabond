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
#include <iostream>

struct vec3;

class Plucker;

/**
 * \class Molecule
 * \brief The highest level group of atoms directly referenced by a Crystal object.
 */

class Molecule : public AtomGroup
{
public:
	Molecule();
	virtual ~Molecule() {};

	virtual void addAtom(AtomPtr atom);
	
	virtual void summary();
	virtual void tieAtomsUp();
	void tiedUpScattering(double *tied, double *all);
	virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal,
	                            int conformer = -1);
	virtual void graph(std::string) {};
	virtual void differenceGraphs(std::string, CrystalPtr) {};
	void setAnchors();

	virtual void refitToSavedPositions();

	virtual void reportParameters();

	double getAbsoluteBFacSubtract()
	{
		return _absoluteBFacSubtract;
	}
	
	static double getAbsoluteBFacSubtract(void *object)
	{
		return static_cast<Molecule *>(object)->_absoluteBFacSubtract;
	}

	void setAbsoluteBFacSubtract(double subtract)
	{
		_absoluteBFacSubtract = subtract;
	}
	
	static void setAbsoluteBFacSubtract(void *object, double subtract);

	double getAbsoluteBFacMult();
	double getAbsoluteBFacSubt();
	void chelate(std::string element, double bufferB = 0);

	void setAbsoluteBFacMult(double mult);

	void setChainID(std::string chain)
	{
		_chainID = chain;
		_name = _chainID;
	}

	std::string getChainID()
	{
		return _chainID;
	}

	virtual std::string getClassName()
	{
		return "Molecule";
	}

	bool isPolymer()
	{
		return (getClassName() == "Polymer");
	}

	bool isWaterNetwork()
	{
		return (getClassName() == "WaterNetwork");
	}

	MoleculePtr shared_from_this()
	{
		AtomGroupPtr groupPtr = AtomGroup::shared_from_this();
		return ToMoleculePtr(groupPtr);
	}
	
	virtual void addParamCounts(int *pos, int *flex)
	{
		
	}
protected:
	virtual std::string getParserIdentifier()
	{
		return "chain_" + _chainID; 
	}

	virtual void addProperties();
	virtual void postParseTidy();    

	void setChangedRotation()
	{
		_changedRotations = true;
	}

private:
	double _absoluteBFacSubtract;
	double _absoluteBFacMult;

	bool _changedRotations;
	std::string _chainID;
};

#endif /* defined(__vagabond__Molecule__) */
