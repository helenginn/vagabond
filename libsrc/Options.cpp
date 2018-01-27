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
double Options::_kick = 0.01;
double Options::_dampen = 0.08;
double Options::_bStart = 1.5;
double Options::_bMult = 0.6;
double Options::_minRes = 0.0;
int Options::_enableTests = false;

Options::Options(int argc, const char **argv)
{
	_manual = false;
	_notify = NULL;
	_globalCount = 0;

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

	_numCycles = 0;
	_tie = true;

	/* Note that argv includes our program name */

	std::cout << std::endl;

	if (argc <= 1)
	{
		std::cout << "Please specify a macromolecular model." << std::endl;
		std::cout << "\te.g., vagabond --with-pdb=xxxx.pdb" << std::endl;
		std::cout << std::endl;
		std::cout << "Alternatively, see all options:" << std::endl;
		std::cout << "\tvagabond --help\n" << std::endl;
        return;
	}

	for (int i = 1; i < argc; i++)
	{
		arguments.push_back(argv[i]);
	}
}

void Options::run()
{
	parse();

    if (arguments.size() <= 1)
    {
        return;
    }

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

	std::cout << "Running in " << (_manual ? "manual" : "automatic") << " mode." << std::endl;

	if (diffractions.size() == 1)
	{
		if (crystals.size() == 1)
		{
			DiffractionPtr data = diffractions[0];

			/* sandbox */
			crystals[0]->writeMillersToFile(data, "pre");
			crystals[0]->getDataInformation(data, 2, 1);

			PolymerPtr chainA = ToPolymerPtr(crystals[0]->molecule("A0"));
			chainA->differenceGraphs("test", crystals[0]);

			if (_tie)
			{
				crystals[0]->setAnchors();
				crystals[0]->tieAtomsUp();
			}

			crystals[0]->tiedUpScattering();

			int count = 0;
			crystals[0]->concludeRefinement(count, data);

			if (!_manual)
			{
				refineAll(RefinementModelPos, 3, &count);
				refineAll(RefinementFlexibility, 20, &count);
				refineAll(RefinementModelPos, 3, &count);
				refineAll(RefinementFine, _numCycles, &count);
			}
			else if (_notify)
			{
				_notify->enable();
			}
		}
	}

	if (!_manual)
	{
		std::cout << std::endl << "**** Finished. ****" << std::endl;
		std::cout << std::endl;
	}
}

void Options::displayHelp()
{
	std::cout << "Syntax: vagabond [options]\n\n" << std::endl;
	std::cout << "Takes an atomistic PDB file and refines it against" << std::endl;
	std::cout << "a reflection list in torsion space.\n\n" << std::endl;
	std::cout << "--help\t\t\t\tDisplays command list.\n" << std::endl;
	std::cout << "--with-pdb=<filename>\t\tName of the input PDB file to refine.\n" << std::endl;
	std::cout << "--with-mtz=<filename>\t\tName of the MTZ file to refine.\n" << std::endl;
	std::cout << "--output-dir=<directory>\tOptional name of a directory to dump processing.\n" << std::endl;
	std::cout << "--anchor-res=<num>\t\tOptional override default anchor residue for all\n\t\t\t\tchains (under development)\n" << std::endl;
	std::cout << "--kick=<num>\t\t\tOptional override for kick fraction for initial bond\n" << std::endl;
	std::cout << "--dampen=<num>\t\t\tOptional override for dampen fraction for all bonds\n" << std::endl;
	std::cout << "--enable-tests=<num>\t\tEnable whatever it is Helen is currently working on\n" << std::endl;
	std::cout << std::endl;
	std::cout << "Baseline command to start running vagabond:\n" << std::endl;
	std::cout << "\tvagabond --with-pdb=start.pdb --with-mtz=start.mtz\n" << std::endl;

	exit(0);
}

void Options::notifyGUI(bool enable)
{
    if (_notify && enable)
    {
        _notify->enable();
    }
    else if (_notify && !enable)
    {
        _notify->disable();
    }
}

void Options::parse()
{
	for (int i = 0; i < arguments.size(); i++)
	{
		bool understood = false;
		std::string arg = arguments[i];

		std::string prefix("--help");

		if (!arg.compare(0, prefix.size(), prefix))
		{
			displayHelp();
		}

		prefix = "--with-pdb=";

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
			diffraction = boost::static_pointer_cast<Diffraction>(mtz);
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

		prefix = "--kick=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string kick_string = arg.substr(prefix.size());
			_kick = atof(kick_string.c_str());

			if (crystals.size() == 0)
			{
				shout_at_user("Initial kick specified, but a coordinate\n"\
							  "file has not been specified yet. Please use\n"\
							  "--with-pdb= to specify some atomic coordinates.");
			}
			else
			{
				CrystalPtr crystal = crystals.at(crystals.size() - 1);
				std::cout << "Setting " << crystal->getFilename()
				<< " to a initial kick of " << _kick << "." << std::endl;

				understood = true;
			}
		}

		prefix = "--bfactor=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string bee_string = arg.substr(prefix.size());
			_bStart = atof(bee_string.c_str());

			understood = true;
		}


		prefix = "--min-res=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string min_string = arg.substr(prefix.size());
			_minRes = atof(min_string.c_str());
                        std::cout << "Minimum resolution set to " << _minRes
                        << " Å." << std::endl;

			understood = true;
		}


		prefix = "--dampen=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string dampen_string = arg.substr(prefix.size());
			_dampen = atof(dampen_string.c_str());

			if (crystals.size() == 0)
			{
				shout_at_user("Overall B factor specified, but a coordinate\n"\
							  "file has not been specified yet. Please use\n"\
							  "--with-pdb= to specify some atomic coordinates.");
			}
			else
			{
				CrystalPtr crystal = crystals.at(crystals.size() - 1);
				std::cout << "Setting " << crystal->getFilename()
				<< " to a dampening of " << _dampen << "." << std::endl;

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

		prefix = "--output-dir=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			_outputDir = arg.substr(prefix.size());
			std::cout << "Setting output directory to ";
			std::cout << _outputDir << "." << std::endl << std::endl;
			FileReader::setOutputDirectory(_outputDir);
			understood = true;
		}
		prefix = "--enable-tests=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string testString = arg.substr(prefix.size());
			_enableTests = atoi(testString.c_str());
			std::cout << "Enabling Helen's test/sandbox no. " <<
			enableTests() << "." << std::endl;
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
			bool madeJoke = parseJoke(arg);

			if (!madeJoke)
			{
				warn_user("I did not understand your command:\n\t"\
						  + arg);
			}
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

void Options::refineAll(RefinementType type, int numCycles, int *count, bool keepGoing)
{
        notifyGUI(false);

	double lastRWork = 200;
	if (count == NULL)
	{
		count = &_globalCount;
	}

	for (int i = 0; i < numCycles; i++)
	{
		for (int j = 0; j < crystals[0]->moleculeCount(); j++)
		{
			MoleculePtr molecule = crystals[0]->molecule(j);
			refinementCycle(molecule, count, type);
		}

		(*count)++;
		double newRWork = crystals[0]->concludeRefinement(*count,
														  diffractions[0]);

		/* Do we go for another cycle? */
		if (keepGoing && i + 1 == numCycles && newRWork < lastRWork)
		{
			std::cout << "Going for another cycle..." << std::endl;
			numCycles++;
		}
		else if (i + 1 == numCycles)
		{
			std::cout << "Leaving it there." << std::endl;
		}

		lastRWork = newRWork;
	}

        notifyGUI(true);
}

void Options::superimposeAll(CrystalPtr crystal)
{
	notifyGUI(false);

	if (!crystal)
	{
		if (crystalCount())
		{
			crystal = crystals[0];
		}
		else return;
	}

	for (int i = 0; i < crystal->moleculeCount(); i++)
	{
		MoleculePtr molecule = crystal->molecule(i);
		if (molecule->isPolymer())
		{
			PolymerPtr polymer = ToPolymerPtr(molecule);
			polymer->superimpose();
			polymer->propagateChange();
			polymer->refreshPositions();
		}
	}

	_globalCount++;
	crystal->concludeRefinement(_globalCount, diffractions[0]);

        notifyGUI(true);
}

void Options::applyBMultiplier()
{
    notifyGUI(false);
    CrystalPtr crystal = crystals[0];

    for (int i = 0; i < crystal->moleculeCount(); i++)
	{
		MoleculePtr molecule = crystal->molecule(i);
		
		if (!molecule->isPolymer())
		{
			std::cout << "Changing B multiplier for HETATMs to: " << _bMult << std::endl;
			molecule->setAbsoluteBFacMult(_bMult);
			molecule->propagateChange();
			molecule->refreshPositions();
		}
	}

	_globalCount++;
	crystal->concludeRefinement(_globalCount, diffractions[0]);

    notifyGUI(true);
}

void Options::refinementCycle(MoleculePtr molecule, int *count,
							  RefinementType type)
{
	DiffractionPtr data = diffractions[0];

	if (molecule->getClassName() == "Polymer")
	{
		PolymerPtr polymer = ToPolymerPtr(molecule);
		polymer->test();

		molecule->refine(crystals[0], type);

		if (_manual) return;

		if (type == RefinementFlexibility)
		{
		//	polymer->optimiseTranslationTensor();
		}

		if (molecule->getClassName() == "Polymer" && (*count == 1))
		{
			polymer->superimpose();
		}
	}
}

bool Options::parseJoke(std::string arg)
{
	if (arg == "vag")
	{
		std::cout << "Unparsed command: vag\n\tI appreciate the offer but let's"\
		"\n\tkeep this vaguely professional\n" << std::endl;
		return true;
	}

	return false;
}

