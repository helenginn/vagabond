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
#include <cstdarg>

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
	void addHydrogens(AtomPtr minor, std::vector<std::string> hNames);
	void addHydrogens(AtomList group, int hNum, ...);
	bool hasHydrogens(BondPtr bond);
	void setNewGeometry(AtomList group, double bondAngle, double torsion,
	                    double portion = -1);
	void setSpin(AtomList group);

	AtomPtr prepareNewHydrogen(AtomPtr parent);
};




#endif
