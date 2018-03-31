//
//  BoneDensity.h
//  vagabond
//
//  Created by Helen Ginn on 31/03/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__BoneDensity__
#define __vagabond__BoneDensity__

#include <map>
#include "shared_ptrs.h"

/**
 * \class BoneDensity
 * \brief Analysis of backbone density and generation of heuristics (maybe)
 * for other parts of the program.
 */

typedef std::map<int, double> DensityScoreMap;

class BoneDensity
{
public:
	BoneDensity();
	void analyse();
	
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;	
	}
	
	void setPolymer(PolymerPtr polymer)
	{
		_polymer = polymer;
	}
private:
	void validate();
	void perMonomerScores();
	void findInflections();
	void createRefinementStrategies();
	DensityScoreMap _densityMap;
	DensityScoreMap _summaryMap;
	
	CrystalPtr _crystal;
	PolymerPtr _polymer;
};


#endif



