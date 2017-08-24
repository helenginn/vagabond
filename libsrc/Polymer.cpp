//
//  Polymer.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Polymer.h"
#include "Sidechain.h"
#include "Monomer.h"
#include "Sampler.h"
#include <iostream>
#include "Backbone.h"
#include "Shouter.h"
#include "Atom.h"
#include "Bond.h"
#include "CSV.h"
#include <fstream>

void Polymer::addMonomer(MonomerPtr monomer)
{
	long existingMonomers = monomerCount();
	_monomers[existingMonomers] = monomer;
	if (monomer)
	{
		monomer->setPolymer(shared_from_this());
	}
}

void Polymer::tieAtomsUp()
{
	for (int i = 0; i < monomerCount(); i++)
	{
		if (getMonomer(i))
		{
			getMonomer(i)->tieAtomsUp();
		}
	}
}

void Polymer::summary()
{
	Molecule::summary();
	std::cout << "| I am a polymer with " << monomerCount() << " monomers." << std::endl;

}

void Polymer::refine(CrystalPtr target, RefinementType rType)
{
	time_t wall_start;
	time(&wall_start);

	for (int i = 0; i < monomerCount(); i++)
	{
		MonomerPtr monomer = getMonomer(i);

		if (!monomer)
		{
			continue;
		}

		BackbonePtr backbone = monomer->getBackbone();

		if (backbone)
		{
			backbone->refine(target, rType);
		}


		SidechainPtr victim = monomer->getSidechain();

		if (victim && victim->canRefine())
		{
			victim->refine(target, rType);
		}
	}

	shout_timer(wall_start, "refinement");

}

void Polymer::makePDB(std::string filename)
{
	std::ofstream file;
	file.open(filename.c_str());

	for (int i = 0; i < monomerCount(); i++)
	{
		MonomerPtr monomer = getMonomer(i);

		if (!monomer)
		{
			continue;
		}

		SidechainPtr victim = monomer->getSidechain();

		file << monomer->getBackbone()->getPDBContribution();
		file << victim->getPDBContribution();
	}

	file.close();

	std::cout << "Written PDB to " << filename << "." << std::endl;
}

void Polymer::graph(std::string graphName)
{
	CSVPtr csv = CSVPtr(new CSV(3, "resnum", "newB", "oldB"));

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		BackbonePtr backbone = getMonomer(i)->getBackbone();
		AtomPtr ca = backbone->findAtom("CA");
		ModelPtr model = ca->getModel();

		if (model->getClassName() != "Bond")
		{
			continue;
		}

		BondPtr bond = std::static_pointer_cast<Bond>(model);
		double meanSq = bond->getMeanSquareDeviation();
		double value = i;

		csv->addEntry(3, value, meanSq, ca->getInitialBFactor());
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = graphName;
	plotMap["height"] = "700";
	plotMap["width"] = "1200";
	plotMap["xHeader0"] = "resnum";
	plotMap["yHeader0"] = "newB";
	plotMap["xHeader1"] = "resnum";
	plotMap["yHeader1"] = "oldB";
	plotMap["colour0"] = "black";
	plotMap["colour1"] = "red";
	plotMap["yMin0"] = "0";
	plotMap["yMin1"] = "0";
	plotMap["yMax0"] = "40";
	plotMap["yMax1"] = "40";

	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "B factor";
	plotMap["style0"] = "line";
	plotMap["style1"] = "line";

	csv->plotPNG(plotMap);
}

