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
#include "DiffractionMTZ.h"
#include "Shouter.h"
#include <iomanip>
#include "Polymer.h"
#include "FileReader.h"

OptionsPtr Options::options;

Options::Options(int argc, const char **argv)
{

	std::cout << "   _______                                _______\n";
	std::cout << " |        ---__________________________---       |\n";
	std::cout << "  \\ o          o   o   o    o   o   o         o /\n";
	std::cout << "    \\ o                                     o /\n";
	std::cout << "      \\ o    \\______           ______/    o /\n";
	std::cout << "        \\ o   \\_____/         \\_____/   o /\n";
	std::cout << "          \\ o           ___           o /\n";
	std::cout << "            \\ o        /   \\        o /\n";
	std::cout << "              -_______-     -_______-\n\n";
	std::cout << "             Vagabond at your service.\n" << std::endl;

	_numCycles = 6;
	_tie = true;

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
		std::cout << "\te.g. --with-mtz=xxxx.mtz\n" << std::endl;
		std::cout << "I will just have a look at your models now." << std::endl;
		std::cout << std::endl;

		outputCrystalInfo();
	}

	if (diffractions.size() == 1)
	{
		if (crystals.size() == 1)
		{
			DiffractionPtr data = diffractions[0];

			/* sandbox */
			crystals[0]->writeCalcMillersToFile(data, "pre");
			crystals[0]->getDataInformation(data, 3, 2);

			if (_tie)
			{
				crystals[0]->setAnchors();
				crystals[0]->tieAtomsUp();
			}

			crystals[0]->tiedUpScattering();

			PolymerPtr chainA = ToPolymerPtr(crystals[0]->molecule("A"));
	//		chainA->downWeightResidues(10, 19, 0.0);
	//		chainA->downWeightResidues(120, 124, 0.0);

			int count = 0;
			crystals[0]->concludeRefinement(count, data, crystals[0]);

			if (_numCycles > 0)
			{
				refineAll(RefinementModelPos, 2, &count);
				refineAll(RefinementFine, _numCycles, &count);
			}
		}
	}

	std::cout << std::endl << "**** Finished. ****" << std::endl;

	std::cout << std::endl;
}

void Options::parse()
{
	for (int i = 0; i < arguments.size(); i++)
	{
		bool understood = false;
		std::string arg = arguments[i];

		std::string prefix("--with-pdb=");

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string pdb_name = arg.substr(prefix.size());

			PDBReader pdb = PDBReader();
			pdb.setFilename(pdb_name);
			CrystalPtr crystal = pdb.getCrystal();

			objects.push_back(crystal);
			crystals.push_back(crystal);
			understood = true;
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
			understood = true;
		}

		prefix = "--target-flex=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string bee_string = arg.substr(prefix.size());
			double bee = atof(bee_string.c_str());

			if (crystals.size() == 0)
			{
				shout_at_user("Overall B factor specified, but a coordinate\n"\
							  "file has not been specified yet. Please use\n"\
							  "--with-pdb= to specify some atomic coordinates.");
			}
			else
			{
				CrystalPtr crystal = crystals.at(crystals.size() - 1);
				crystal->setOverallBFactor(bee);
				std::cout << "Setting " << crystal->getFilename()
				<< " to a target flexiness of " << bee << "." << std::endl;
				understood = true;
			}
		}

		prefix = "--anchor-res=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			if (crystals.size() == 0)
			{
				shout_at_user("Anchor residue specified, but a coordinate\n"\
							  "file has not been specified yet. Please use\n"\
							  "--with-pdb= to specify some atomic coordinates.");
			}

			std::string anchor_string = arg.substr(prefix.size());

			size_t comma = anchor_string.find(",");
			CrystalPtr crystal = crystals.at(crystals.size() - 1);

			while (true)
			{
				std::string number = anchor_string;

				if (comma < anchor_string.size())
				{
					number = anchor_string.substr(0, comma);
					anchor_string = anchor_string.substr(comma + 1);
					comma = anchor_string.find(",");
				}

				int anchor = atoi(number.c_str());
				crystal->addAnchorResidue(anchor);

				if (comma >= anchor_string.size())
				{
					int anchor = atoi(anchor_string.c_str());
					crystal->addAnchorResidue(anchor);
					break;
				}
			}

			understood = true;
		}

		prefix = "--max-res=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string max_res = arg.substr(prefix.size());
			double maxRes = atof(max_res.c_str());

			if (crystals.size() == 0)
			{
				shout_at_user("Max resolution specified, but a coordinate\n"\
							  "file has not been specified yet. Please use\n"\
							  "--with-pdb= to specify some atomic coordinates.");
			}
			else
			{
				CrystalPtr crystal = crystals.at(crystals.size() - 1);
				crystal->setMaxResolution(maxRes);
				std::cout << "Setting " << crystal->getFilename()
				<< " to resolution " << maxRes << " Å." << std::endl;
				understood = true;
			}
		}

		prefix = "--num-cycles=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string result = arg.substr(prefix.size());
			_numCycles = atoi(result.c_str());
			understood = true;
		}

		prefix = "--no-tie";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			_tie = false;
			understood = true;
		}

		if (!understood)
		{
			warn_user("I did not understand your command:\n\t"\
					  + arg);
		}
	}
}

void Options::outputCrystalInfo()
{
	if (!crystals.size())
	{
		return;
	}

	shout_at_user("Sorry I do actually need a data set.");

	for (int i = 0; i < crystals.size(); i++)
	{
		crystals[i]->realSpaceClutter();
		crystals[i]->fourierTransform(1);
	}
}

void Options::refineAll(RefinementType type, int numCycles, int *count)
{
	for (int i = 0; i < numCycles; i++)
	{
		for (int j = 0; j < crystals[0]->moleculeCount(); j++)
		{
			MoleculePtr molecule = crystals[0]->molecule(j);
			refinementCycle(molecule, count, type);
		}

		(*count)++;
		crystals[0]->concludeRefinement(*count, diffractions[0], crystals[0]);
	}

}

void Options::refinementCycle(MoleculePtr molecule, int *count,
							  RefinementType type)
{
	DiffractionPtr data = diffractions[0];

	if (molecule->getClassName() == "Polymer")
	{
		PolymerPtr polymer = ToPolymerPtr(molecule);
		polymer->test();

		if (true)
		{
			molecule->refine(crystals[0], type);
		}

		if (molecule->getClassName() == "Polymer" && (*count == 1))
		{
			polymer->superimpose();
		}

		if (molecule->getClassName() == "Polymer" && (*count == 2))
		{
	//		crystals[0]->changeAnchors(1);
		}
	}

	/*
	 else if (molecule->getClassName() == "Polymer" && i == 2)
	 {
	 if (crystals[0]->totalAnchors() <= 1) continue;

	 count++;

	 polymer->scaleFlexibilityToBFactor(crystals[0]);
	 crystals[0]->concludeRefinement(count, data, crystals[0]);
	 }*/
	/*					else if (molecule->getClassName() == "Polymer" && i % 2 == 0)
	 {
	 int myAnchor = (i / 2) - 1;
	 if (myAnchor >= crystals[0]->totalAnchors()) continue;

	 count++;
	 crystals[0]->changeAnchors(myAnchor);
	 crystals[0]->concludeRefinement(count, data, crystals[0]);
	 }*/
}
