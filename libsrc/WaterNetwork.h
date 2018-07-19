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

/**
 * \class WaterNetwork
 * \brief A subclass of Molecule which looks after all HOH heteroatoms,
 * and any other water molecules generated.
 */

class WaterNetwork : public Molecule
{
public:
	WaterNetwork();
	virtual ~WaterNetwork() {}
	
	virtual void summary();

	virtual std::string getClassName()
	{
		return "WaterNetwork";
	}
	
	virtual void addProperties();
	
	void partitionNetworks(CrystalPtr crystal);
	
	void setMonomer(MonomerPtr _monomer);
	static double vsRefineWaterNetwork(void *object);
private:	
	std::vector<WaterClusterPtr> _clusters;
	void reportOnClusters();
};

#endif
