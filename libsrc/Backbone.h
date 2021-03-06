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

/**
 * \class Backbone
 * \brief A class containing references to the atoms in a backbone of a
 * single Monomer, including the shared atom between the Monomer and the
 * Sidechain.
 */



class Backbone : public AtomGroup
{
public:
	virtual ~Backbone() {};
	
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
		return (_timesRefined > 2);
	}

	AtomPtr betaCarbonTorsionAtom();

	virtual void refine(CrystalPtr target, RefinementType rType);
	void setAnchor();
protected:
	virtual bool shouldRefineAtom(AtomPtr atom);
	virtual void addProperties();

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
