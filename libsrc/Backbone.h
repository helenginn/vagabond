//
//  Backbone.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Backbone__
#define __vagabond__Backbone__

#include <stdio.h>
#include "AtomGroup.h"
#include "Sampler.h"

class Backbone : public AtomGroup
{
public:
	void setPolymer(PolymerPtr poly)
	{
		_myPolymer = poly;
	}
	
	PolymerPtr getPolymer()
	{
		return _myPolymer.lock();
	}

	void setResNum(int resNum)
	{
		_resNum = resNum;
	}

	int getResNum()
	{
		return _resNum;
	}

	AtomPtr betaCarbonTorsionAtom();

	virtual void refine(CrystalPtr target, RefinementType rType);
	void setAnchor();
protected:
	virtual bool shouldRefineMagicAxis(BondPtr bond);
private:
	int _resNum;
	PolymerWkr _myPolymer;
};

#endif /* defined(__vagabond__Backbone__) */
