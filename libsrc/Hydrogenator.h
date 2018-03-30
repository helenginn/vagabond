//
// Hydrogenator.h
// vagabond
//
// Created by Helen Ginn on 22/03/2018
// Copyright (c) 2018 Helen Ginn
//

#ifndef __vagabond__Hydrogenator__
#define __vagabond__Hydrogenator__

#include "shared_ptrs.h"

class Hydrogenator
{
public:
	Hydrogenator();
	
	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}
	
	void hydrogenate();

private:
	MonomerPtr _monomer;
};




#endif
