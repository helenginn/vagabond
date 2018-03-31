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
#include "VBondReader.h"
#include "SSRigger.h"

OptionsPtr Options::options;
double Options::_kick = 0.01;
int Options::_solvent = 1;
double Options::_dampen = 0.08;
double Options::_bStart = 1.5;
double Options::_bMult = 1;
double Options::_minRes = 0.0;
int Options::_enableTests = 3;
bool Options::_powder = false;
double Options::_sampling = -1;

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

	if (!_manual && !crystals.size())
	{
		shout_at_user("I need a model.\n"\
		              "\te.g. --with-pdb=xxxx.pdb");
	}

	if (!_manual && !datasets.size())
	{
		std::cout << "Warning: You have not specified any data sources!" << std::endl;
		std::cout << "\te.g. --with-mtz=xxxx.mtz\n" << std::endl;
		std::cout << "I will just have a look at your models now." << std::endl;
		std::cout << std::endl;

		outputCrystalInfo();
	}

	std::cout << std::setprecision(3) << std::endl;
	std::cout << "Running in " << (_manual ? "manual" : "automatic") << " mode." << std::endl;

	if (diffractions.size() == 1)
	{
		if (crystals.size() == 1)
		{
			DiffractionPtr data = diffractions[0];

			CrystalPtr crystal = getActiveCrystal();

			/* sandbox */
			crystal->writeMillersToFile(data, "pre");
			crystal->getDataInformation(data, 2, 1);

			if (_tie)
			{
				crystal->setAnchors();
				crystal->tieAtomsUp();
			}

			crystal->tiedUpScattering();

			int count = 0;
			crystal->concludeRefinement(count, data);

			if (shouldPowder())
			{
				crystal->makePowders();
				goto finished;
			}

			if (!_manual)
			{
				refineAll(RefinementModelPos, 3, &count);
				refineAll(RefinementModelPos, 3, &count);
				refineAll(RefinementFine, _numCycles, &count);
			}
			else if (_notify)
			{
				_notify->enable();
			}
		}
	}

	finished:

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

		prefix = "--with-vbond=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string vbond_name = arg.substr(prefix.size());

			VBondReader vReader = VBondReader();
			vReader.setFilename(vbond_name);
			CrystalPtr crystal = vReader.getCrystal();
			//            PolymerPtr polymer = ToPolymerPtr(crystal->molecule(0));
			//            polymer->splitConformers();

			understood = true;

			if (!crystal)
			{
				std::cout << "Read failed." << std::endl;
			}
			else
			{
				objects.push_back(crystal);
				crystals.push_back(crystal);
			}
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

		prefix = "--solvent=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string solvent_string = arg.substr(prefix.size());
			_solvent = atoi(solvent_string.c_str());
			understood = true;
		}

		prefix = "--kick=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string kick_string = arg.substr(prefix.size());
			_kick = atof(kick_string.c_str());
			understood = true;
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

		prefix = "--powder";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			_powder = true;
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

		statusMessage("Calculating R factors...");

		double newRWork = getActiveCrystal()->concludeRefinement(*count,
		                                                         diffractions[0]);

		agreementSummary();

		/* Do we go for another cycle? */
		if (keepGoing && i + 1 == numCycles && newRWork < lastRWork)
		{
			std::cout << "Going for another cycle..." << std::endl;
			numCycles++;
		}

		lastRWork = newRWork;
	}

	notifyGUI(true);
}

void Options::superimposeAll(CrystalPtr crystal)
{
	notifyGUI(false);

	statusMessage("Superimposing complete ensemble...");

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
			polymer->refreshPositions(true);
		}
	}

	statusMessage("Calculating R factors...");

	recalculateFFT();

	notifyGUI(true);
}

void Options::applyBMultiplier()
{
	notifyGUI(false);

	CrystalPtr crystal = getActiveCrystal();

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

	recalculateFFT();

	notifyGUI(true);
}

void Options::findDisulphides()
{
	notifyGUI(false);

	CrystalPtr crystal = getActiveCrystal();

	statusMessage("Finding disulphide bonds.");

	SSRigger rigger;
	rigger.setCrystal(crystal);
	rigger.findDisulphides();
	
	recalculateFFT();

	notifyGUI(true);
}

void Options::openModel(std::string pdbName)
{
	if (crystals.size())
	{
		statusMessage("Already have a crystal, not loading another. Not yet.");
		return;
	}

	ModelFile modelType = ModelFilePDB;

	// if the name ends in 'vbond', let's switch...
	if (pdbName.substr(pdbName.length() - 5, 5) == "vbond")
	{
		modelType = ModelFileVagabond;
	}

	notifyGUI(false);

	if (modelType == ModelFilePDB)
	{
		statusMessage("Loading PDB file " + pdbName + "...");
	}
	else if (modelType == ModelFileVagabond)
	{
		statusMessage("Loading Vagabond file " + pdbName + "...");
	}

	CrystalPtr crystal = CrystalPtr();

	if (modelType == ModelFilePDB)
	{ 
		PDBReader pdb = PDBReader();
		pdb.setFilename(pdbName);
		crystal = pdb.getCrystal();
	}
	else if (modelType == ModelFileVagabond)
	{
		VBondReader reader = VBondReader();
		reader.setFilename(pdbName);
		crystal = reader.getCrystal();
	}

	objects.push_back(crystal);
	crystals.push_back(crystal);

	if (modelType == ModelFilePDB && _tie)
	{
		statusMessage("Tying up atoms...");
		crystals[0]->setAnchors();
		crystals[0]->tieAtomsUp();
	}

	for (int i = 0; i < crystal->moleculeCount(); i++)
	{
		MoleculePtr molecule = crystal->molecule(i);
		if (molecule->isPolymer())
		{
			PolymerPtr polymer = ToPolymerPtr(molecule);
			polymer->refreshPositions(true);
		}
	}

	crystals[0]->tiedUpScattering();

	if (diffractions.size())
	{
		recalculateFFT();
	}
	else if (modelType == ModelFilePDB)
	{
		statusMessage("Loaded PDB file " + pdbName + ".");
	}
	else if (modelType == ModelFileVagabond)
	{
		statusMessage("Loaded Vagabond file " + pdbName + ".");
	}

	notifyGUI(true);
}

void Options::openMTZ(std::string mtzName)
{
	if (diffractions.size())
	{
		statusMessage("Sorry, diffraction data has already loaded. Please restart!");
		return;
	}

	notifyGUI(false);

	statusMessage("Loading MTZ file " + mtzName + "...");

	DiffractionMtzPtr mtz;
	mtz = DiffractionMtzPtr(new DiffractionMtz());
	DiffractionPtr diffraction;
	diffraction = boost::static_pointer_cast<Diffraction>(mtz);
	diffraction->setFilename(mtzName);
	diffraction->load();

	objects.push_back(diffraction);
	datasets.push_back(diffraction);
	diffractions.push_back(diffraction);

	if (crystals.size())
	{
		recalculateFFT();
	}
	else
	{
		statusMessage("Loaded MTZ file " + mtzName + ".");
	}

	notifyGUI(true);
}

std::string Options::rTypeString(RefinementType type)
{
	switch (type)
	{
		case RefinementModelPos:
		return "Torsions against PDB positions";
		case RefinementFine:
		return "Torsions against electron density";
		default:
		return "Unknown";
	}
}

void Options::refinementCycle(MoleculePtr molecule, int *count,
                              RefinementType type)
{
	DiffractionPtr data = diffractions[0];

	std::string rString = rTypeString(type);
	statusMessage("Refining structure, method: " + rString);

	if (molecule->getClassName() == "Polymer")
	{
		PolymerPtr polymer = ToPolymerPtr(molecule);
		polymer->test();

		statusMessage("Refining structure, chain " + molecule->getChainID() + "...");
		molecule->refine(crystals[0], type);

		if (_manual) return;

		if (molecule->getClassName() == "Polymer" && (*count == 1))
		{
			polymer->superimpose();
		}
	}
}

void Options::fitWholeMolecule(bool translation, bool rotation)
{
	notifyGUI(false);

	statusMessage("Applying whole-molecule fits...");

	CrystalPtr crystal = getActiveCrystal();
	crystal->fitWholeMolecules(translation, rotation);

	statusMessage("Ready.");
	notifyGUI(true);
}

void Options::recalculateFFT()
{
	if (!diffractions.size()) return;
	if (!crystals.size()) return;

	_globalCount++;
	DiffractionPtr mtz = diffractions[0];
	statusMessage("Calculating R factors...");
	crystals[0]->concludeRefinement(_globalCount, mtz);
	agreementSummary();
}

void Options::statusMessage(std::string message)
{
	if (_notify)
	{
		_notify->setMessage(message);
	}

	std::cout << message << std::endl;
}

void Options::agreementSummary()
{
	notifyGUI(false);    

	statusMessage(crystals[0]->agreementSummary());

	notifyGUI(true);
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


