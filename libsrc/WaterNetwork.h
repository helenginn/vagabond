//
//  WaterNetwork.h
//  vagabond
//
//  Created by Helen Ginn on 01/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__WaterNetwork__
#define __vagabond__WaterNetwork__

#include <stdio.h>
#include "shared_ptrs.h"
#include "Molecule.h"
#include <vector>
#include <map>
#include "Options.h"

typedef struct
{
	vec3 *left;
	vec3 *middle;
	double weight;
	double ratio; /* left-to-middle over middle-to-right */

} WaterCalc;

typedef std::vector<WaterCalc> WaterCalcs;

/**
 * \class WaterNetwork
 * \brief A subclass of Molecule which looks after all HOH heteroatoms,
 * and any other water molecules generated.
 * 
 */

class WaterNetwork : public Molecule
{
public:
	WaterNetworkPtr shared_from_this()
	{
		return ToWaterNetworkPtr(Parser::shared_from_this());
	}

	WaterNetwork();
	virtual ~WaterNetwork() {}
	
	virtual void summary();

	virtual std::string getClassName()
	{
		return "WaterNetwork";
	}
	
	void recalculate();
	void calculateSingle(SpongePtr sponge, bool others = false);
	void refreshSponges();
	
	double score();
	static double sScore(void *object)
	{
		return static_cast<WaterNetwork *>(object)->score();
	}
	
	virtual void addProperties();
	virtual void postParseTidy();
	
	void setMonomer(MonomerPtr _monomer);
private:	
	void setActive(int a);
	void generateCalculations();
	SpongePtr findFirstSponge();
	void refineSponges();
	std::vector<SpongePtr> acquireSponges(SpongePtr seed);

	std::vector<SpongePtr> _sponges;
	std::vector<WaterCalcs> _calcs;

	int _n;
};

#endif

