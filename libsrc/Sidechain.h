//
//  Sidechain.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Sidechain__
#define __vagabond__Sidechain__

#include <stdio.h>
#include "AtomGroup.h"
#include "Sampler.h"
#include "FileReader.h"

/**
 * \class Sidechain
 * \brief A class containing references to the atoms in a sidechain of a
 * single Monomer, including the shared atom between the Monomer and the
 * Backbone.
 */


class Sidechain : public AtomGroup
{
public:
	Sidechain()
	{
		_resNum = 0;
		_canRefine = true;
		_rotamerised = false;
	}
	
	virtual ~Sidechain() {};

	bool canRefine()
	{
		return _canRefine;
	}

	void setCanRefine(bool canRefine)
	{
		_canRefine = canRefine;
	}

	void setResNum(int resNum)
	{
		_resNum = resNum;
	}

	void setPolymer(PolymerPtr poly)
	{
		_myPolymer = poly;
	}

	PolymerPtr getPolymer()
	{
		return _myPolymer.lock();
	}

	virtual bool shouldRefineAngles()
	{
		return (_timesRefined > 2);
	}
	
	bool isRotamerised()
	{
		return _rotamerised;	
	}

	void splitConformers(int count);
	void parameteriseAsRotamers();
	virtual void refine(CrystalPtr target, RefinementType rType);


protected:
	virtual AtomList topLevelAtoms()
	{
		return findAtoms("CB");
	}

	virtual std::string getClassName()
	{
		return "Sidechain";
	}

	virtual std::string getParserIdentifier()
	{
		return "side_" + i_to_str(_resNum);
	}

	virtual void addProperties();
private:
	bool _rotamerised;
	bool _canRefine;
	int _resNum;
	PolymerWkr _myPolymer;
};

#endif /* defined(__vagabond__Sidechain__) */
