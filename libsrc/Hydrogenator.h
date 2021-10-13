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

/**
 * \class Hydrogenator
 * \brief Adds hydrogens to a model and sets the geometry to default values
 * for each relevant atom, if hydrogens are not already present. */

class Hydrogenator
{
public:
	Hydrogenator();
	
	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}
	
	void hydrogenate();
	void addHydrogens(AtomList group, int hNum, ...);

	/** Pass negative torsion if you do not want it set */
	void setNewGeometry(AtomList group, double bondAngle, 
	                    double torsion = -1,
	                    double portion = -1);

	static void adjustProlineHydrogens(BondPtr newBond, bool current);
private:

	void getHydrogenBonds(AtomGroupPtr group);
	void addHydrogens(AtomPtr minor, std::vector<std::string> hNames);
	bool hasHydrogens(BondPtr bond);
	void setSpin(AtomList group);
	void adjustBond(BondPtr newBond);
	AtomPtr prepareNewHydrogen(AtomPtr parent);
	
	double getHBondLength(AtomPtr minor);
	
	double _cAlpha2H;
	MonomerPtr _monomer;

};




#endif
