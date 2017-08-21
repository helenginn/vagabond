//
//  AtomGroup.h
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__AtomGroup__
#define __vagabond__AtomGroup__

#include <stdio.h>
#include <string>
#include <vector>
#include "shared_ptrs.h"

class AtomGroup
{
public:
	AtomPtr findAtom(std::string atomType);
	AtomList findAtoms(std::string atomType);

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

	void addBond(BondPtr bond)
	{
		_bonds.push_back(bond);
	}

	int bondCount()
	{
		return _bonds.size();
	}

	BondPtr bond(int i)
	{
		return _bonds[i];
	}

	double totalElectrons();

	std::string getPDBContribution();

	void setTied()
	{
		_beenTied = true;
	}
protected:
	AtomGroup();
	int _timesRefined;

	bool isTied()
	{
		return _beenTied;
	}
private:
	MonomerWkr _monomer;

	std::vector<BondPtr> _bonds;
	std::vector<AtomPtr> _atoms;

	bool _beenTied;
private:
	
};

#endif /* defined(__vagabond__AtomGroup__) */
