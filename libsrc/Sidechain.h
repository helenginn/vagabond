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

class Sidechain : public AtomGroup, public Sampler
{
public:
	Sidechain()
	{
		_canRefine = false;
	}

	void refine(CrystalPtr target, RefinementType rType);
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

	void fixBackboneTorsions(AtomPtr betaTorsion);

private:
	bool _canRefine;
	int _resNum;
	PolymerWkr _myPolymer;
};

#endif /* defined(__vagabond__Sidechain__) */
