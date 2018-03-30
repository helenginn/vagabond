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
#include "FileReader.h"

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

	virtual bool shouldRefineAngles()
	{
		return (_timesRefined > 0);
	}

	AtomPtr betaCarbonTorsionAtom();

	virtual void refine(CrystalPtr target, RefinementType rType);
	void setAnchor();
protected:
	virtual void addProperties();
	virtual bool shouldRefineMagicAxis(BondPtr);

	virtual std::string getClassName()
	{
		return "Backbone";
	}

	virtual std::string getParserIdentifier()
	{
		return "back_" + i_to_str(_resNum);
	}

private:
	int _resNum;
	PolymerWkr _myPolymer;
};

#endif /* defined(__vagabond__Backbone__) */
