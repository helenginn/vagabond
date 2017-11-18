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
#include "Sampler.h"
#include <map>

class AtomGroup : public boost::enable_shared_from_this<AtomGroup>, public Sampler
{
public:
	AtomPtr findAtom(std::string atomType);
	AtomPtr findAtom(std::string atomType, std::string confID);
	AtomList findAtoms(std::string atomType);

	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}

	MonomerPtr getMonomer()
	{
		return _monomer.lock();
	}

	virtual void addAtom(AtomPtr atom)
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
		return _bonds[i].lock();
	}

	double totalElectrons();
	double getAverageBFactor(bool initial = false);
	double getAverageDisplacement();

	std::string getPDBContribution(PDBType pdbType,
								   CrystalPtr crystal = CrystalPtr());

	void setTied()
	{
		_beenTied = true;
	}

	void setUseAbsolute();
	int totalElectrons(int *fcWeighted);

	virtual void refine(CrystalPtr target, RefinementType rType);
	void setWeighting(double value);
	void resetMagicAxes();
	int conformerCount();
	std::string conformer(int i);
	void propagateChange();
protected:
	AtomGroup();
	void addAtomsFrom(AtomGroupPtr child);
	virtual AtomList topLevelAtoms();
	bool hasAtom(AtomPtr anAtom);

	bool isTied()
	{
		return _beenTied;
	}
private:
	MonomerWkr _monomer;

	std::vector<BondWkr> _bonds;
	std::vector<AtomPtr> _atoms;

	bool _beenTied;

	std::map<std::string, int> conformerMap();

};

#endif /* defined(__vagabond__AtomGroup__) */
