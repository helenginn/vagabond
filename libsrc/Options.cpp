//
//  Options.cpp
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Options.h"
#include "VScript.h"
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
double Options::_kick = 0.005;
int Options::_solvent = 1;
double Options::_dampen = 0.08;
double Options::_bStart = 1.5;
double Options::_bMult = 1;
double Options::_minRes = -1.0;
double Options::_maxRes = -1.0;
int Options::_enableTests = 3;
bool Options::_powder = false;
double Options::_sampling = -1;
std::string Options::_solventFile;

Options::Options(int argc, const char **argv)
{
	_parsed = false;
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
	if (!_parsed)
	{
		parse();
	}
	
	if (_minRes < 0) _minRes = 0;

	if (_mtzFile.length())
	{
		DiffractionMtzPtr mtz;
		mtz = DiffractionMtzPtr(new DiffractionMtz());
		DiffractionPtr diffraction;
		diffraction = boost::static_pointer_cast<Diffraction>(mtz);
		diffraction->setFilename(_mtzFile);
		diffraction->load();

		objects.push_back(diffraction);
		datasets.push_back(diffraction);
		diffractions.push_back(diffraction);
	}
	
	int pdbFile = -1;

	if (_modelFile.length())
	{
		if (_modelFile.substr(_modelFile.length() - 3, 3) == "pdb")
		{
			std::cout << "Guessing PDB file" << std::endl;
			pdbFile = 1;
		}	
		else
		{
			std::cout << "Guessing vbond file" << std::endl;
			pdbFile = 0;
		}
	}

	if (pdbFile == 1)
	{
		PDBReader pdb = PDBReader();
		pdb.setFilename(_modelFile);
		CrystalPtr crystal = pdb.getCrystal();

		if (!crystal)
		{
			std::cout << "Read of " << _modelFile << " failed." << std::endl;
		}
		else
		{
			objects.push_back(crystal);
			crystals.push_back(crystal);
		}
	}
	else
	{
		VBondReader vReader = VBondReader();
		vReader.setFilename(_modelFile);
		CrystalPtr crystal = vReader.getCrystal();

		objects.push_back(crystal);
		crystals.push_back(crystal);
	}
	
	if (_outputDir.length())
	{
		FileReader::setOutputDirectory(_outputDir);
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
			if (_tie)
			{
				crystal->setAnchors();
				crystal->tieAtomsUp();
			}

			crystal->tiedUpScattering();

			recalculateFFT();
			
			if (_notify)
			{
				_notify->setInstruction(InstructionTypeResetExplorer);
			}

			if (shouldPowder())
			{
				crystal->makePowders();
				recalculateFFT();
				goto finished;
			}
			
			executeScript();

			if (_notify)
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
	
	FFT::cleanupPlans();
}

void Options::executeScript()
{
	if (!_scriptName.length())
	{
		return;
	}

	VScript script = VScript();
	std::string contents = get_file_contents(_scriptName);
	
	script.loadScript(contents);
	script.execute();
}

void Options::displayHelp()
{
	std::cout << "Syntax: vagabond [options]\n\n" << std::endl;
	std::cout << "Takes an atomistic PDB file and refines it against" << std::endl;
	std::cout << "a reflection list in torsion space.\n\n" << std::endl;
	std::cout << "--help\t\t\t\tDisplays command list.\n" << std::endl;
	std::cout << "--with-model=<filename>\t\tName of the input PDB or Vagabond model file to refine.\n" << std::endl;
	std::cout << "--with-mtz=<filename>\t\tName of the MTZ file to refine against.\n" << std::endl;
	std::cout << "--output-dir=<directory>\tOptional name of a directory to dump processing.\n" << std::endl;
	std::cout << "--kick=<num>\t\t\tOptional override for kick fraction for initial bond\n" << std::endl;
	std::cout << "--dampen=<num>\t\t\tOptional override for dampen fraction for all bonds\n" << std::endl;
	std::cout << "--bfactor=<num>\t\t\tOptional override for the assigned B factor for anchor residues\n" << std::endl;
	std::cout << "--min-res=<value>\t\tOptional override for the minimum resolution in Ångströms.\n" << std::endl;
	std::cout << "--max-res=<value>\t\tOptional override for the maximum resolution in Ångströms.\n" << std::endl;
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

int Options::parseParameter(std::string arg, std::string prefix, double *ptr)
{
	if (!arg.compare(0, prefix.size(), prefix))
	{
		std::string val_string = arg.substr(prefix.size());
		*ptr = atof(val_string.c_str());
		std::cout << "Parsing " << prefix << ", setting to " << *ptr
		<< " units." << std::endl;

		return true;
	}
	
	return false;
}

int Options::parseParameter(std::string arg, std::string prefix,
                            std::string *ptr)
{
	if (!arg.compare(0, prefix.size(), prefix))
	{
		std::string val_string = arg.substr(prefix.size());
		*ptr = val_string;
		std::cout << "Parsing " << prefix << ", setting to " << *ptr
		<< std::endl;

		return true;
	}
	
	return false;
}

void Options::parse()
{
	for (size_t i = 0; i < arguments.size(); i++)
	{
		int understood = false;
		std::string arg = arguments[i];

		std::string prefix("--help");

		if (!arg.compare(0, prefix.size(), prefix))
		{
			displayHelp();
		}

		understood |= parseParameter(arg, "--with-vscript=", &_scriptName);
		understood |= parseParameter(arg, "--solvent-file=", &_solventFile);

		understood |= parseParameter(arg, "--with-model=", &_modelFile);
		understood |= parseParameter(arg, "--with-mtz=", &_mtzFile);

		prefix = "--solvent=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string solvent_string = arg.substr(prefix.size());
			_solvent = atoi(solvent_string.c_str());
			understood = true;
		}

		understood |= parseParameter(arg, "--min-res=", &_minRes);
		understood |= parseParameter(arg, "--bfactor=", &_bStart);
		understood |= parseParameter(arg, "--kick=", &_kick);
		understood |= parseParameter(arg, "--dampen=", &_dampen);
		understood |= parseParameter(arg, "--max-res=", &_maxRes);
		understood |= parseParameter(arg, "--output-dir=", &_outputDir);

		prefix = "--enable-tests=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string testString = arg.substr(prefix.size());
			_enableTests = atoi(testString.c_str());
			std::cout << "Enabling Helen's test/sandbox no. " <<
			enableTests() << "." << std::endl;
			understood = true;
		}

		prefix = "--solvent-analysis";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			_powder = true;
			understood = true;
		}

		prefix = "--diff-matrix=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string diff = arg.substr(prefix.size());
			_diffMatrix = diff;
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
	
	_parsed = true;
}

void Options::outputCrystalInfo()
{
	if (!crystals.size())
	{
		return;
	}

	shout_at_user("Sorry I do actually need a data set.");
}

void Options::backboneAnalysis()
{
	notifyGUI(false);

	getActiveCrystal()->backboneDensityAnalysis();
	statusMessage("Refining via backbone-derived plan.");

	notifyGUI(true);
	statusMessage("Ready.");
}

void Options::refineAll(RefinementType type, int numCycles, int *count, bool)
{
	notifyGUI(false);

	if (count == NULL)
	{
		count = &_globalCount;
	}

	for (int i = 0; i < numCycles; i++)
	{
		for (size_t j = 0; j < crystals[0]->moleculeCount(); j++)
		{
			MoleculePtr molecule = crystals[0]->molecule(j);
			refinementCycle(molecule, type);
		}

		statusMessage("Calculating R factors...");

		recalculateFFT();
	}

	notifyGUI(true);
}

void Options::diffMatrix()
{
	std::string contents = get_file_contents(_diffMatrix);
	std::vector<std::string> lines = split(contents, '\n');
	
	std::vector<vec3> pos1s, pos2s;
	
	for (int i = 0; i < lines.size(); i++)
	{
		std::vector<std::string> components = split(lines[i], ',');
		
		if (components.size() < 6) continue;
		
		vec3 a; vec3 b;
		a.x = atof(components[0].c_str());
		a.y = atof(components[1].c_str());
		a.z = atof(components[2].c_str());
		b.x = atof(components[3].c_str());
		b.y = atof(components[4].c_str());
		b.z = atof(components[5].c_str());
		
		pos1s.push_back(a);
		pos2s.push_back(b);
	}
	
	std::cout << "i, j, value" << std::endl;
	
	for (int i = 0; i < pos1s.size(); i++)
	{
		for (int j = 0; j < pos1s.size(); j++)
		{
			vec3 a = pos1s[i];
			vec3 other_a = pos1s[j];
			vec3 b = pos2s[j];
			vec3 other_b = pos2s[i];

			vec3 diff1 = vec3_subtract_vec3(a, b);
			vec3 diff2 = vec3_subtract_vec3(other_a, other_b);
			
			double length1 = vec3_length(diff1);	
			double length2 = vec3_length(diff2);	

			std::cout << i << ", " << j << ", " << length2 - length1;
			std::cout << std::endl;
		}
		
	}

	std::cout << std::endl;
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

	for (size_t i = 0; i < crystal->moleculeCount(); i++)
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

	for (size_t i = 0; i < crystal->moleculeCount(); i++)
	{
		MoleculePtr molecule = crystal->molecule(i);

		if (!molecule->isPolymer())
		{
			std::cout << "Changing B multiplier for HETATMs to: " << _bMult << std::endl;
			molecule->setAbsoluteBFacMult(_bMult);
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
		getActiveCrystal()->setAnchors();
		getActiveCrystal()->tieAtomsUp();
	}

	getActiveCrystal()->hydrogenateContents();

	for (size_t i = 0; i < crystal->moleculeCount(); i++)
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
		case RefinementSidechain:
		return "Sidechains only against electron density";
		default:
		return "Unknown";
	}
}

void Options::refinementCycle(MoleculePtr molecule, RefinementType type)
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

void Options::recalculateFFT(bool saveState)
{
	if (!diffractions.size()) return;
	if (!crystals.size()) return;

	_globalCount++;
	DiffractionPtr mtz = diffractions[0];
	statusMessage("Calculating R factors...");
	
	getActiveCrystal()->concludeRefinement(_globalCount, mtz);
	agreementSummary();
	
	if (saveState)
	{
		getActiveCrystal()->saveState();
	}

	std::cout << "Total states: " << getActiveCrystal()->stateCount() <<
	std::endl;
}

void Options::previousState()
{
	int state = -1;
	
	if (crystals[0]->stateCount() == 1)
	{
		state = 0;	
	}
	
	statusMessage("Restoring previous state...");
	crystals[0]->restoreState(state);
	
	recalculateFFT(false);
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


double Options::getActiveCrystalDStar()
{
	CrystalPtr crystal = getActiveCrystal();
	DiffractionPtr data = getActiveData();
	
	return crystal->getMaximumDStar(data);
}
