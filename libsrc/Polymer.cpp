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
	CSVPtr csv = CSVPtr(new CSV(2, "resnum", "rmsd"));

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

		csv->addEntry(2, value, meanSq);
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = graphName;
	plotMap["height"] = "700";
	plotMap["width"] = "1000";
	plotMap["xHeader0"] = "resnum";
	plotMap["yHeader0"] = "rmsd";

	plotMap["xTitle0"] = "Residue number";
	plotMap["yTitle0"] = "RMSD";
	plotMap["style0"] = "line";

	csv->plotPNG(plotMap);
}

