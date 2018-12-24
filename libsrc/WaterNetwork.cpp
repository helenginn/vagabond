//
//  WaterNetwork.cpp
//  vagabond
//
//  Created by Helen Ginn on 01/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "WaterNetwork.h"
#include "WaterCluster.h"
#include "Crystal.h"
#include "Atom.h"

#define MAX_HBOND_DISTANCE 3.4

WaterNetwork::WaterNetwork()
{

}

void WaterNetwork::summary()
{
	Molecule::summary();
	std::cout << "| I am a water network with " << atomCount() 
	<< " waters." << std::endl;
}

void WaterNetwork::addProperties()
{
	Molecule::addProperties();
}

bool hasWaters(AtomGroupPtr list)
{
	for (int i = 0; i < list->atomCount(); i++)
	{
		if (list->atom(i)->isHeteroAtom() &&
		    list->atom(i)->getAtomName() == "O")
		{
			return true;
		}
	}
	
	return false;
}

