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
#include "Monomer.h"

class Sidechain
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

	long atomCount()
	{
		return _atoms.size();
	}

	AtomPtr atom(int i)
	{
		return _atoms[i];
	}

	AtomPtr findAtom(std::string atomType);
private:
	MonomerWkr _monomer;

	std::vector<AtomPtr> _atoms;
};

#endif /* defined(__vagabond__Sidechain__) */
