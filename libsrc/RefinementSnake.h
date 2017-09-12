//
//  RefinementSnake.h
//  vagabond
//
//  Created by Helen Ginn on 11/09/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__RefinementSnake__
#define __vagabond__RefinementSnake__

#include <stdio.h>
#include "RefinementStrategy.h"

class RefinementSnake : public RefinementStrategy
{
public:
	RefinementSnake() : RefinementStrategy()
	{
		
	};

	void addBond(BondPtr bond)
	{
		_bonds.push_back(bond);
	}

	long bondCount()
	{
		return _bonds.size();
	}

	BondPtr bond(int i)
	{
		return _bonds[i];
	}

	void setParentSampler(Sampler *sampler);
	virtual void refine();

private:
	std::vector<BondPtr> _bonds;
	std::vector<AtomPtr> _atoms;
	Sampler *_parentSampler;

	std::vector<double> getTorsionList();
	void applyTorsionList(std::vector<double> list);
};

#endif /* defined(__vagabond__RefinementSnake__) */
