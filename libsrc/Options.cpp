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
#include "Polymer.h"

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
	std::cout << "Vagabond at your service.\n" << std::endl;

	parse();

	if (!crystals.size())
	{
		shout_at_user("I need a model.\n"\
					  "\te.g. --with-pdb=xxxx.pdb");
	}

	if (!datasets.size())
	{
		std::cout << "Warning: You have not specified any data sources!" << std::endl;
		std::cout << "\te.g. --with-mtz=xxxx.mtz\n" << std::endl;
		std::cout << "I will just have a look at your models now." << std::endl;
		std::cout << std::endl;

		outputCrystalInfo();
	}

	if (diffractions.size() == 1)
	{
		if (crystals.size() == 1)
		{
			int prop = 3;
			/* sandbox */
			DiffractionPtr data = diffractions[0];
			crystals[0]->realSpaceClutter();
			crystals[0]->transplantAmplitudes(data, prop, prop-1);
			MoleculePtr molecule = crystals[0]->molecule("A");
			/*
			for (int i = 0; i < 1; i++)
			{
				molecule->refine(crystals[0], RefinementBroad);
				crystals[0]->realSpaceClutter();
				crystals[0]->transplantAmplitudes(data, prop, prop-1);
			}
*/
			for (int i = 0; i < 4; i++)
			{
				molecule->refine(crystals[0], RefinementFine);
				crystals[0]->realSpaceClutter();
				crystals[0]->transplantAmplitudes(data, prop, prop-1);
			}

			crystals[0]->fourierTransform(1);
			crystals[0]->writeCalcMillersToFile();

		}
	}

	std::cout << std::endl;
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
		crystals[i]->realSpaceClutter();
		crystals[i]->fourierTransform(1);
		crystals[i]->writeCalcMillersToFile();
	}
}