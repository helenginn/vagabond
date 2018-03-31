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
#include <vector>
#include "shared_ptrs.h"
#include "Sampler.h"

typedef struct
{
	int startRes;
	int endRes;
	RefinementType rType;
} BackboneInstruction;

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
	
	size_t instructionCount()
	{
		return _instructions.size();
	}
	
	BackboneInstruction instruction(int i)
	{
		return _instructions[i];
	}
private:
	void validate();
	void perMonomerScores();
	void findInflections();
	void createRefinementStrategies();
	DensityScoreMap _densityMap;
	DensityScoreMap _summaryMap;
	
	std::vector<BackboneInstruction> _instructions;
	
	CrystalPtr _crystal;
	PolymerPtr _polymer;
};


#endif



