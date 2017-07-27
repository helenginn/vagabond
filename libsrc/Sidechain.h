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

	void refine(CrystalPtr target);
	bool canRefine()
	{
		return _canRefine;
	}

	void setCanRefine(bool canRefine)
	{
		_canRefine = canRefine;
	}
private:
	bool _canRefine;
};

#endif /* defined(__vagabond__Sidechain__) */
