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
#include "Monomer.h"

class Backbone
{
public:
	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}

	MonomerPtr getMonomer()
	{
		return _monomer.lock();
	}

	void addAtom(AtomPtr atom)
	{
		_atoms.push_back(atom);
	}
	
private:
	MonomerWkr _monomer;

	std::vector<AtomPtr> _atoms;
};

#endif /* defined(__vagabond__Backbone__) */
