//
//  Options.cpp
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Options.h"
#include <iostream>
#include "PDBReader.h"
#include "Crystal.h"
#include "DiffractionMtz.h"
#include "Shouter.h"
#include <iomanip>

Options::Options(int argc, const char **argv)
{
	/* Note that argv includes our program name */

	std::cout << std::endl;

	if (argc <= 1)
	{
		std::cout << "Please specify a macromolecular model." << std::endl;
		std::cout << "\te.g., vagabond --with-pdb=xxxx.pdb" << std::endl;
		std::cout << std::endl;
		exit(0);
	}

	for (int i = 1; i < argc; i++)
	{
		arguments.push_back(argv[i]);
	}
}

void Options::run()
{
	parse();

	if (!crystals.size())
	{
		shout_at_user("I need a model.\n"\
					  "\te.g. --with-pdb=xxxx.pdb");
	}

	if (!datasets.size())
	{
		std::cout << "Warning: You have not specified any data sources!" << std::endl;
		std::cout << "         At the moment, you can't anyway." << std::endl;
		std::cout << std::endl;

		outputCrystalInfo();
	}

	if (diffractions.size() == 1)
	{
		if (crystals.size() == 1)
		{
			crystals[0]->scaleToDiffraction(diffractions[0]);
			double rFac = crystals[0]->rFactorWithDiffraction(diffractions[0]);

			std::cout << "Crystal has R factor with data of " <<
			std::setprecision(3) << rFac * 100 << "%." << std::endl;
			crystals[0]->transplantAmplitudes(diffractions[0]);
			crystals[0]->writeCalcMillersToFile();
		}
	}
}

void Options::parse()
{
	for (int i = 0; i < arguments.size(); i++)
	{
		std::string arg = arguments[i];

		std::string prefix("--with-pdb=");

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string pdb_name = arg.substr(prefix.size());

			PDBReader pdb = PDBReader();
			pdb.setFilename(pdb_name);
			ObjectPtr crystal = std::static_pointer_cast<Object>(pdb.getCrystal());

			objects.push_back(crystal);
			crystals.push_back(pdb.getCrystal());
		}

		prefix = "--with-mtz=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string mtz_name = arg.substr(prefix.size());

			DiffractionMtzPtr mtz;
			mtz = DiffractionMtzPtr(new DiffractionMtz());
			DiffractionPtr diffraction;
			diffraction = std::static_pointer_cast<Diffraction>(mtz);
			diffraction->setFilename(mtz_name);
			diffraction->load();

			objects.push_back(diffraction);
			datasets.push_back(diffraction);
			diffractions.push_back(diffraction);
		}
	}
}

void Options::outputCrystalInfo()
{
	if (!crystals.size())
	{
		return;
	}

	std::cout << "Converting " << crystals.size() <<
	" pdbs to reflection list." << std::endl << std::endl;

	for (int i = 0; i < crystals.size(); i++)
	{
		crystals[i]->calculateMillers();
		crystals[i]->writeCalcMillersToFile();
	}
}