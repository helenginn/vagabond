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
	CSVPtr csvDamp = CSVPtr(new CSV(4, "resnum", "dN-CA", "dCA-C", "dC-N"));
	CSVPtr csvBlur = CSVPtr(new CSV(4, "resnum", "bN-CA", "bCA-C", "bC-N"));

	for (int i = 0; i < monomerCount(); i++)
	{
		if (!getMonomer(i))
		{
			continue;
		}

		BackbonePtr backbone = getMonomer(i)->getBackbone();
		AtomPtr ca = backbone->findAtom("CA");
		ModelPtr caModel = ca->getModel();
		AtomPtr n = backbone->findAtom("N");
		ModelPtr nModel = n->getModel();
		AtomPtr c = backbone->findAtom("C");
		ModelPtr cModel = c->getModel();

		double value = i;
		double caDampen = 0; double cDampen = 0; double nDampen = 0;
		double caBlur = 0; double cBlur = 0; double nBlur = 0;

		if (caModel->getClassName() == "Bond")
		{
			BondPtr caBond = std::static_pointer_cast<Bond>(caModel);
			double meanSq = caBond->getMeanSquareDeviation();
			csv->addEntry(3, value, meanSq, ca->getInitialBFactor());
			caDampen = Bond::getDampening(&*caBond);
			caBlur = Bond::getTorsionBlur(&*caBond);

			if (caDampen > 0) caBlur = 0;
		}
		else
		{
			csv->addEntry(3, value, ca->getInitialBFactor(), ca->getInitialBFactor());
		}

		if (cModel->getClassName() == "Bond")
		{
			BondPtr cBond = std::static_pointer_cast<Bond>(cModel);
			cDampen = Bond::getDampening(&*cBond);
			cBlur = Bond::getTorsionBlur(&*cBond);

			if (cDampen > 0) cBlur = 0;
		}

		if (nModel->getClassName() == "Bond")
		{
			BondPtr nBond = std::static_pointer_cast<Bond>(nModel);
			nDampen = Bond::getDampening(&*nBond);
			nBlur = Bond::getTorsionBlur(&*nBond);
			if (nDampen > 0) nBlur = 0;
		}

		csvDamp->addEntry(4, value, caDampen, cDampen, nDampen);
		csvBlur->addEntry(4, value, caBlur, cBlur, nBlur);
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

	plotMap["filename"] = "dampening_" + graphName;
	plotMap["yHeader0"] = "dN-CA";
	plotMap["xHeader1"] = "resnum";
	plotMap["yHeader1"] = "dCA-C";
	plotMap["xHeader2"] = "resnum";
	plotMap["yHeader2"] = "dC-N";
	plotMap["style2"] = "line";
	plotMap["colour2"] = "blue";
	plotMap["yTitle0"] = "Dampening factor";
	plotMap["yMin0"] = "-2";
	plotMap["yMin1"] = "-2";
	plotMap["yMin2"] = "-2";
	plotMap["yMax0"] = "2";
	plotMap["yMax1"] = "2";
	plotMap["yMax2"] = "2";

	csvDamp->plotPNG(plotMap);

	plotMap["filename"] = "blurring_" + graphName;
	plotMap["yHeader0"] = "bN-CA";
	plotMap["yHeader1"] = "bCA-C";
	plotMap["yHeader2"] = "bC-N";
	plotMap["yTitle0"] = "Blurring factor";
	plotMap["yMin0"] = "-0.5";
	plotMap["yMin1"] = "-0.5";
	plotMap["yMin2"] = "-0.5";
	plotMap["yMax0"] = "0.5";
	plotMap["yMax1"] = "0.5";
	plotMap["yMax2"] = "0.5";

	csvBlur->plotPNG(plotMap);
}

