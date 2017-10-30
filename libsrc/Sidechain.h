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

class Sidechain : public AtomGroup
{
public:
	Sidechain()
	{
		_canRefine = false;
	}

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

	void setInitialDampening();
	void fixBackboneTorsions(AtomPtr betaTorsion);
protected:
	virtual bool shouldRefineMagicAxis(BondPtr bond);
	virtual AtomPtr topLevelAtom()
	{
		return findAtom("CB");
	}
private:
	bool _canRefine;
	int _resNum;
	PolymerWkr _myPolymer;
};

#endif /* defined(__vagabond__Sidechain__) */
