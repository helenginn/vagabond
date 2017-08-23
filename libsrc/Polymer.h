//
//  Polymer.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Polymer__
#define __vagabond__Polymer__

#include <stdio.h>
#include "shared_ptrs.h"
#include "Molecule.h"
#include <vector>
#include <map>

class Polymer : public Molecule, public std::enable_shared_from_this<Polymer>
{
public:
	void addMonomer(MonomerPtr monomer);
	virtual void summary();
	virtual void tieAtomsUp();
	virtual void refine(CrystalPtr target, RefinementType rType);
	virtual void makePDB(std::string filename);
	void graph(std::string graphName);

	void addUnknownMonomers(int number)
	{
		for (int i = 0; i < number; i++)
		{
			addMonomer(MonomerPtr());
		}
	}

	MonomerPtr getMonomer(int i)
	{
		return _monomers[i];
	}

	long monomerCount()
	{
		return _monomers.size();
	}

	std::string className()
	{
		return "Polymer";
	}
private:
	std::map<long, MonomerPtr> _monomers;

};

#endif /* defined(__vagabond__Polymer__) */
