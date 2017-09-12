//
//  Polymer.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Polymer__
#define __vagabond__Polymer__

#include <stdio.h>
#include "shared_ptrs.h"
#include "Molecule.h"
#include "Sampler.h"
#include <vector>
#include <map>

class Polymer :
public Molecule,
public std::enable_shared_from_this<Polymer>,
public Sampler
{
public:
	Polymer()
	{
		_dampening = 0.05;
		_anchorNum = 0;
	}

	void addMonomer(MonomerPtr monomer);
	virtual void summary();
	virtual void tieAtomsUp();
	virtual void refine(CrystalPtr target, RefinementType rType);
	virtual void makePDB(std::string filename);
	void graph(std::string graphName);

	static double getConstantDampening(void *object);
	static void setConstantDampening(void *object, double value);

	static void setInitialKick(void *object, double value);
	static double getInitialKick(void *object);

	void scaleFlexibilityToBFactor(double value);

	void changeAnchor(int num);
	void setAnchor(int num)
	{
		_anchorNum = num;
	}

	int getAnchor()
	{
		return _anchorNum;
	}

	void addUnknownMonomers(int number)
	{
		for (int i = 0; i < number; i++)
		{
			addMonomer(MonomerPtr());
		}
	}

	MonomerPtr getMonomer(int i)
	{
		return _monomers[i];
	}

	long monomerCount()
	{
		return _monomers.size();
	}

	virtual std::string getClassName()
	{
		return "Polymer";
	}
private:
	void refineMonomer(MonomerPtr monomer, CrystalPtr target,
					   RefinementType rType);

	std::map<long, MonomerPtr> _monomers;

	int _anchorNum;
	double _dampening;
};

#endif /* defined(__vagabond__Polymer__) */
