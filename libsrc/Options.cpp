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
#include "DiffractionMTZ.h"
#include "Shouter.h"
#include <iomanip>
#include "Polymer.h"
#include "FileReader.h"

Options::Options(int argc, const char **argv)
{
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
			double propFo = 4;
			double propFc = 3.2;

			/* sandbox */
			DiffractionPtr data = diffractions[0];

			MoleculePtr molecule = crystals[0]->molecule("A");
			PolymerPtr polymer = ToPolymerPtr(molecule);
/*
			polymer->downWeightResidues(856, 861, 0);
			polymer->downWeightResidues(883, 890, 0);
			polymer->downWeightResidues(1009, 1015, 0);
*/
			crystals[0]->writeCalcMillersToFile(data, "pre");
			crystals[0]->getDataInformation(data, propFo, propFc);

			if (_tie)
			{
				crystals[0]->setAnchors();
				crystals[0]->tieAtomsUp();
			}

			crystals[0]->tiedUpScattering();


			int count = 0;
			crystals[0]->concludeRefinement(count, data, crystals[0]);

			RefinementType type = RefinementModelRMSD;

			if (_numCycles > 0)
			{
				for (int i = 0; i < 25; i++)
				{
					if (true)
					{
						count++;
						molecule->refine(crystals[0], type);
						crystals[0]->concludeRefinement(count, data, crystals[0]);
					}

					if (molecule->getClassName() == "Polymer" && i == 0)
					{
                        count++;
						polymer->minimiseRotations();
						polymer->minimiseCentroids();
						crystals[0]->concludeRefinement(count, data, crystals[0]);
					}
/*
					else if (molecule->getClassName() == "Polymer" && i == 2)
					{
						if (crystals[0]->totalAnchors() <= 1) continue;

						count++;
						crystals[0]->changeAnchors(1);
						polymer->scaleFlexibilityToBFactor(crystals[0]);
						crystals[0]->concludeRefinement(count, data, crystals[0]);
					}*/
					else if (molecule->getClassName() == "Polymer" && i == 8)
					{
						count++;
						polymer->minimiseRotations();
						polymer->minimiseCentroids();
						crystals[0]->concludeRefinement(count, data, crystals[0]);
					}
/*					else if (molecule->getClassName() == "Polymer" && i % 2 == 0)
					{
						int myAnchor = (i / 2) - 1;
						if (myAnchor >= crystals[0]->totalAnchors()) continue;

						count++;
						crystals[0]->changeAnchors(myAnchor);
						crystals[0]->concludeRefinement(count, data, crystals[0]);
                    }*/

				}
			}

			for (int i = 0; i < _numCycles; i++)
			{
				count++;
				if (i >= 10)
				{
					molecule->refine(crystals[0], RefinementFineBlur);
				}
				else
				{
					molecule->refine(crystals[0], RefinementFine);
				}

				crystals[0]->concludeRefinement(count, data, crystals[0]);
			}
		}
	}

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

