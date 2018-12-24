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
#include "SSRigger.h"
#include "DiffractionMTZ.h"
#include "Shouter.h"
#include <iomanip>
#include "Polymer.h"
#include "WaterNetwork.h"
#include "FileReader.h"
#include "VBondReader.h"

OptionsPtr Options::options;
double Options::_kick = 0.000;
int Options::_solvent = 1;
int Options::_nCycles = 5;
double Options::_dampen = 0.0;
double Options::_bStart = 20.;
double Options::_bMult = 1;
double Options::_bSubt = -1;
double Options::_bReal = -1;
double Options::_minRes = 0.0;
double Options::_maxRes = -1.0;
double Options::_probeRadius = -0.1;
bool Options::_useRFree = true;
int Options::_enableTests = 3;
int Options::_bondAngles = 1;
bool Options::_powder = false;
bool Options::_overfit = false;
bool Options::_diagnostics = false;

bool Options::_refine = false;
bool Options::_far = false;
bool Options::_rPosition = true;
bool Options::_rSidechains = true;
bool Options::_rInter = true;
bool Options::_rIntra = true;
bool Options::_peptideMovement = true;

std::string Options::_anchor = "";
ScalingType Options::_scaleType = ScalingTypeShell;
double Options::_sampling = -1;
int Options::_nSamples = -1;
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
		std::cout << "\te.g., vagabond --with-model=xxxx.pdb" << std::endl;
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
	
	writeCommandLine();
	
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
	
	openModel(_modelFile);
	
	notifyGUI(false);

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

	std::cout << std::setprecision(3);

	if (crystals.size() >= 1)
	{
		CrystalPtr crystal = getActiveCrystal();

		if (_tie)
		{
			crystal->setAnchors();
			crystal->tieAtomsUp();
			getActiveCrystal()->hydrogenateContents();
		}

		crystal->tiedUpScattering();

		if (_notify)
		{
			_notify->setInstruction(InstructionTypeResetExplorer);
		}
		
		if (_overfit)
		{
			crystal->overfitTest();
		}

		executeProtocol();

		executeScript();

		if (_notify)
		{
			_notify->enable();
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

void Options::executeProtocol()
{
	if (!_refine)
	{
		std::cout << "No refinement protocol selected." << std::endl;
		recalculateFFT();
		return;
	}
	
	CrystalPtr crystal = getActiveCrystal();
	
	for (int i = 0; i < 5 && _rPosition; i++)
	{
		std::cout << "Refining positions to PDB (" << 
		i + 1 << " / 5)" << std::endl;
		crystal->refinePositions();
		
		if (i == 1)
		{
			SSRigger rigger;
			rigger.setCrystal(crystal);
			rigger.findDisulphides();
		}
	}
	
	for (int i = 0; i < 3 && _far; i++)
	{
		recalculateFFT();
		std::cout << "Refining positions to density (" << 
		i + 1 << " / 3)" << std::endl;
		crystal->refineCrude();
	}
	
	recalculateFFT();
	
	for (int i = 0; i < _nCycles; i++)
	{
		if (_rInter || _rIntra)
		{
			std::cout << "Flex macrocycle (" << 
			i + 1 << " / " << _nCycles<< ")" << std::endl;
		}

		if (_rInter)
		{
			double oldWork = crystal->getRWork();
			crystal->fitWholeMolecules();
			recalculateFFT();
			
			double newWork = crystal->getRWork();

			if (newWork > oldWork && !_rIntra)
			{
				crystal->restoreState(-1);
			}
		}

		/* In case we need to do remedial work */
		double oldWork = crystal->getRWork();
		double oldB = crystal->averageBFactor();
		
		for (int i = 0; i < 3 && _rIntra; i++)
		{
			std::cout << "Intramolecular flex microcycle (" << 
			i + 1 << " / 3)" << std::endl;
			crystal->refineIntraMovements();
			recalculateFFT();
		}
		
		double newWork = crystal->getRWork();
		
		if (newWork > oldWork)
		{
			std::cout << "Remedial work to reduce overall flexibility." 
			<< std::endl;
			/* Remedial action required. */
			double newB = crystal->averageBFactor();
			double ratio = oldB / newB;
			/* Only do a smaller increment at a time */
			double increment = (ratio - 1) / 2 + 1;

			bool undone = false;

			while (!undone)
			{
				crystal->scaleAnchorBs(ratio);
				recalculateFFT();
				double latest = crystal->getRWork();
				
				if (latest > newWork)
				{
					std::cout << "Undoing last change and resuming next cycle." 
					<< std::endl;
					crystal->restoreState(-1);
					break;
				}
				else
				{
					newWork = latest;
				}
			}
		}
		else
		{
			std::cout << "Rwork has reduced due to intramolecular flex." 
			<< std::endl;

		}
	}
	
	if (_rSidechains)
	{
		std::cout << "Refining sidechains to density (" << 
		1 << " / 1)" << std::endl;
		getActiveCrystal()->refineSidechains();
		recalculateFFT();
	}
}

void Options::executeScript()
{
	if (!_scriptName.length())
	{
		return;
	}

	VScript script = VScript();
	std::string contents;
	
	try
	{
		contents = get_file_contents(_scriptName);
	}
	catch (int e)
	{
		warn_user("Could not load vscript file: " + _scriptName);
		return;
	}
	
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
//	std::cout << "--with-vscript=<filename>\tName of Vagabond script to run after loading data/model.\n" << std::endl;
	std::cout << "--output-dir=<directory>\tOptional name of a directory to dump processing.\n" << std::endl;
	std::cout << "--kick=<num>\t\t\tOptional override for kick fraction for initial bond\n" << std::endl;
	std::cout << "--dampen=<num>\t\t\tOptional override for dampen fraction for all bonds\n" << std::endl;
	std::cout << "--bfactor=<num>\t\t\tOptional override for the assigned B factor for anchor residues\n" << std::endl;
	std::cout << "--min-res=<value>\t\tOptional override for the minimum resolution in Ångströms.\n" << std::endl;
	std::cout << "--max-res=<value>\t\tOptional override for the maximum resolution in Ångströms.\n" << std::endl;
	std::cout << "--global-b=<value>\t\tAdd a real space B factor when adding explicit atoms to the map.\n"\
	"\t\t\t\t (default 0.).\n" << std::endl;
	std::cout << "--shell-scale\t\t\tWhen calculating R factors, scale each resolution bin of Fcalc"\
	"to Fobs\n\t\t\t\t (default 0).\n" << std::endl;
	std::cout << std::endl;
	std::cout << "--no-rfree\t\t\tDo not use this function.\n" << std::endl;
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

int Options::parseParameter(std::string arg, std::string prefix, int *ptr)
{
	if (!arg.compare(0, prefix.size(), prefix))
	{
		std::string val_string = arg.substr(prefix.size());
		*ptr = atoi(val_string.c_str());
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

int Options::parseParameter(std::string arg, std::string prefix,
                            bool *ptr)
{
	std::string pos = "--" + prefix;
	std::string neg = "--no-" + prefix;

	if (!arg.compare(0, pos.size(), pos))
	{
		std::string val_string = arg.substr(pos.size());
		*ptr = true;
		std::cout << "Setting property " + prefix + " to Yes."
		<< std::endl;

		return true;
	}
	else if (!arg.compare(0, neg.size(), neg))
	{
		std::string val_string = arg.substr(neg.size());
		*ptr = false;
		std::cout << "Setting property " + prefix + " to No."
		<< std::endl;

		return true;
	}
	
	return false;
	
}

void Options::writeCommandLine()
{
	std::cout << "Used command line: " << std::endl << std::endl;
	std::cout << "vagabond";
	
	if (_notify)
	{
		std::cout << "-gui";
	}
	
	std::cout << " ";

	for (int i = 0; i < arguments.size(); i++)
	{
		std::cout << arguments[i] << " ";
	}

	std::cout << std::endl;	
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
		understood |= parseParameter(arg, "--with-script=", &_scriptName);

		understood |= parseParameter(arg, "--with-model=", &_modelFile);
		understood |= parseParameter(arg, "--with-mtz=", &_mtzFile);

		prefix = "--solvent=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string solvent_string = arg.substr(prefix.size());
			_solvent = atoi(solvent_string.c_str());
			understood = true;
		}

		understood |= parseParameter(arg, "--max-res=", &_maxRes);
		understood |= parseParameter(arg, "--min-res=", &_minRes);
		understood |= parseParameter(arg, "--bfactor=", &_bStart);
		understood |= parseParameter(arg, "--bond-angles=", &_bondAngles);
		understood |= parseParameter(arg, "--nsamples=", &_nSamples);
		understood |= parseParameter(arg, "--probe-radius=", &_probeRadius);
		understood |= parseParameter(arg, "--ncycles=", &_nCycles);
		understood |= parseParameter(arg, "--global-b=", &_bReal);
		understood |= parseParameter(arg, "--kick=", &_kick);
		understood |= parseParameter(arg, "--dampen=", &_dampen);
		understood |= parseParameter(arg, "--output-dir=", &_outputDir);
		understood |= parseParameter(arg, "--anchor=", &_anchor);

		prefix = "--enable-tests=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string testString = arg.substr(prefix.size());
			_enableTests = atoi(testString.c_str());
			std::cout << "Enabling Helen's test/sandbox no. " <<
			enableTests() << "." << std::endl;
			understood = true;
		}

		prefix = "--overfit-test";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			_overfit = true;
			understood = true;
		}

		prefix = "--diff-matrix=";

		if (!arg.compare(0, prefix.size(), prefix))
		{
			std::string diff = arg.substr(prefix.size());
			_diffMatrix = diff;
			understood = true;
		}

		understood |= parseParameter(arg, "tie", &_tie);
		understood |= parseParameter(arg, "refine", &_refine);
		understood |= parseParameter(arg, "far", &_far);
		understood |= parseParameter(arg, "position", &_rPosition);
		understood |= parseParameter(arg, "inter-mol", &_rInter);
		understood |= parseParameter(arg, "intra-mol", &_rIntra);
		understood |= parseParameter(arg, "sidechain", &_rSidechains);
		understood |= parseParameter(arg, "rfree", &_useRFree);
		understood |= parseParameter(arg, "diagnostics", &_diagnostics);
		
		int shellNum = 0;
		bool shellstood = parseParameter(arg, "--shell-scale=", &shellNum);
		
		understood |= shellstood;
		
		if (shellstood)
		{
			_scaleType = ScalingType(shellNum);
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

void Options::applyBMultiplier()
{
	notifyGUI(false);

	CrystalPtr crystal = getActiveCrystal();
	
	if (!crystal)
	{
		return;
	}

	std::cout << "Applied HETATM B multiplier " << _bMult <<
	" and subtraction " << _bSubt << "." << std::endl;

	for (size_t i = 0; i < crystal->moleculeCount(); i++)
	{
		MoleculePtr molecule = crystal->molecule(i);

		if (!molecule->isPolymer())
		{
			molecule->forceModelRecalculation();
		}
	}
	

	notifyGUI(true);
}

void Options::findDisulphides()
{
	notifyGUI(false);

	CrystalPtr crystal = getActiveCrystal();

	statusMessage("Finding disulphide bonds.");

	
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
	
	crystal->summary();

	objects.push_back(crystal);

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

	}

	statusMessage("Refining structure, chain " + molecule->getChainID() + "...");
	molecule->refine(crystals[0], type);
}

void Options::recalculateFFT(bool saveState)
{
	statusMessage("Calculating R factors...");
	
	Crystal::vsConcludeRefinement(&*(getActiveCrystal()));
	agreementSummary();
	
//	std::cout << "Total states: " << getActiveCrystal()->stateCount() <<
//	std::endl;
}

void Options::openInCoot()
{
	CrystalPtr crystal = getActiveCrystal();
	crystal->openInCoot();
}

void Options::previousState()
{
	int state = -1;
	
	statusMessage("Undo: restoring previous state...");
	crystals[0]->restoreState(state);
}

void Options::statusMessage(std::string message, bool std_out)
{
	OptionsPtr opt = getRuntimeOptions();

	if (opt->_notify)
	{
		opt->_notify->setMessage(message);
	}

	if (std_out)
	{
		std::cout << message << std::endl;
	}
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

void Options::renderDensity()
{
	if (!_notify)
	{
		return;
	}
	
	_notify->setRenderDensity();
}

void Options::adjustBFactor()
{
	CrystalPtr crystal = getActiveCrystal();
	double adjusted = crystal->getAdjustBFactor();
	double newB = crystal->getRealBFactor();
	
	statusMessage("Adjusted real B factor by " + f_to_str(adjusted, 2)
	              + " to " + f_to_str(newB, 2) + ".");
}

void Options::focusOnPosition(vec3 pos)
{
	if (_notify)
	{
		_notify->focusOnPosition(pos);
	}
}

void Options::pauseGUIFishing(bool on)
{
	if (getRuntimeOptions()->_notify)
	{
		getRuntimeOptions()->_notify->pause(on);
	}
}

void Options::chelate()
{
	getActiveCrystal()->chelate();
}

void Options::omitScan()
{
	getActiveCrystal()->omitScan();
}

void Options::refitBackbone(int start, int end)
{
	MoleculePtr molecule = getActiveCrystal()->molecule(2);
	
	if (!molecule || !molecule->isPolymer())
	{
		molecule = getActiveCrystal()->molecule(0);
	}
	
	PolymerPtr poly = ToPolymerPtr(molecule);
	poly->refitBackbone(start, end);
}
