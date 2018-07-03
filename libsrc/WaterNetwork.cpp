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

void WaterNetwork::refine(CrystalPtr crystal, RefinementType type)
{
	partitionNetworks(crystal);
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

void WaterNetwork::reportOnClusters()
{
	std::cout << "Partitioned into " << _clusters.size() << 
	" hydrogen-bonding clusters." << std::endl;
	
	double total = 0;
	double waters = 0;

	for (int i = 0; i < _clusters.size(); i++)
	{
		int count = _clusters[i]->atomCount();
		int water = 0;
		
		for (int j = 0; j < count; j++)
		{
			AtomPtr atom = _clusters[i]->atom(j);
			
			if (atom->isHeteroAtom() && atom->getAtomName() == "O")
			{
				water++;
			}
		}
		
		total += count;
		waters += water;
		WaterClusterPtr cluster = _clusters[i];
		
		{
			std::cout << " Cluster " << i << " has ";
			std::cout << count << " atoms of which ";
			std::cout << water << " are water.";
			
			if (water < 4)
			{
				std::cout << "This water is no. " <<
				_clusters[i]->findAtoms("O")[0].lock()->getAtomNum() << std::endl;
			}
			
			std::cout << std::endl;
		}
	}
	
	total /= double(_clusters.size());
	waters /= double(_clusters.size());
	
	std::cout << "Each cluster has an average of " << waters << " waters out";
	std::cout << " of " << total << " hydrogen bonders." << std::endl;
	
	std::cout << std::endl << "Now investigating one cluster ..." << std::endl;
	
	_clusters[27]->findNeighbours();
}

void WaterNetwork::partitionNetworks(CrystalPtr crystal)
{
	std::vector<WaterClusterPtr> groups;
	std::vector<AtomPtr> remaining = getAtoms();
	std::vector<AtomPtr> others = crystal->getHydrogenBonders();
	
	std::cout << "Partitioning " << remaining.size() << " waters, ";
	std::cout << "using " << others.size() << " other hydrogen bonders";
	std::cout << " from the crystal." << std::endl;
	
	remaining.reserve(remaining.size() + others.size());
	remaining.insert(remaining.end(), others.begin(), others.end());
	
	while (remaining.size() > 0)
	{
		WaterClusterPtr group = WaterClusterPtr (new WaterCluster());
		group->addAtom(remaining[0]);
		remaining.erase(remaining.begin());

		for (int j = 0; j < group->atomCount(); j++)
		{
			AtomPtr existing = group->atom(j);

			if (!existing->isHeteroAtom())
			{
				continue;
			}

			for (int i = 0; i < remaining.size(); i++)
			{
				AtomPtr poss = remaining[i];

				if (poss->getDistanceFrom(&*existing)
				    < MAX_HBOND_DISTANCE)
				{
					group->addAtom(remaining[i]);
					remaining.erase(remaining.begin() + i);
					i = 0;
				}
			}
		}

		if (hasWaters(group))	
		{
			groups.push_back(group);
		}
	}
	
	_clusters = groups;
	reportOnClusters();
	
}

