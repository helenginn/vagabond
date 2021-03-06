//
//  PDBReader.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "PDBReader.h"
#include "FileReader.h"
#include "vec3.h"
#include "shared_ptrs.h"
#include <vector>
#include "Shouter.h"
#include "mat3x3.h"
#include "Crystal.h"
#include "Polymer.h"
#include "WaterNetwork.h"
#include "Monomer.h"
#include "Absolute.h"
#include <sstream>
#include <iomanip>
#include "Element.h"
#include <unistd.h>

PDBReader::PDBReader()
{
	_rWork = -1;
	_rFree = -1;
	_foundCrystal = false;
	_breakMe = false;
}

AbsolutePtr PDBReader::makeAbsolute(std::string line)
{
	if (line.length() < 78)
	{
		shout_at_user("This ATOM line in the PDB\n"\
		              "is less than 78 characters and might be\n"\
		              "non-standard. Please check the formatting\n"\
		              "of this line.\n\n" + line);
	}

	std::string xData, yData, zData, element, bFactor, occupancy;
	std::string resNum, chainID, atomName, atomNum, resName, alternate;

	atomNum = line.substr(6, 5);
	atomName = line.substr(11, 5);
	
	std::string trName = atomName;
	trim(trName);
	if (_type.length() && trName != _type)
	{
		return AbsolutePtr();
	}

	alternate = line.substr(16, 1);
	resName = line.substr(17, 3);
	chainID = line.substr(21, 1);
	resNum = line.substr(22, 7);
	xData = line.substr(30, 8);
	yData = line.substr(38, 8);
	zData = line.substr(46, 8);
	occupancy = line.substr(54, 6);
	bFactor = line.substr(60, 6);
	element = line.substr(76, 2);

	double xValue = atof(xData.c_str());
	double yValue = atof(yData.c_str());
	double zValue = atof(zData.c_str());
	double bFacValue = atof(bFactor.c_str());
	double occValue = atof(occupancy.c_str());
	int resNumValue = atoi(resNum.c_str());
	double atomNumValue = atoi(atomNum.c_str());

	vec3 vec = make_vec3(xValue, yValue, zValue);
	AbsolutePtr abs = AbsolutePtr(new Absolute(vec, bFacValue, element, occValue));
	_myAbsolute = abs;
	abs->setIdentity(resNumValue, chainID, resName, atomName, atomNumValue);

	if (alternate[0] != ' ')
	{
		abs->setAlternativeConformerName(alternate);
	}

	if (line.substr(0, 6) == "HETATM")
	{
		abs->setHeteroAtom(true);
		checkNotModifiedHetatm(abs);
	}

	return abs;
}

void PDBReader::checkNotModifiedHetatm(AbsolutePtr abs)
{
	if (abs->getResName() == "mly" || abs->getResName() == "sch" ||
	    abs->getResName() == "mlz" || abs->getResName() == "mse")
	{
		Shouter *shout;
		shout = new Shouter("PDB file contains modified (non-standard)\n"\
		                    "amino acid residue as HETATM entry.\n"\
		                    "Please relabel this residue (no. " +
		                    i_to_str(abs->getResNum()) + ") as\n"\
		                    "an ATOM record.");
		throw shout;
	}
}

/* We want to check we should be appending to the correct molecule */
void PDBReader::validateMolecule(AbsolutePtr atom)
{
	if (atom->getChainID() != _myChain || _breakMe)
	{
		/* We have switched proteins. Have we gone back to an old one? */

		/* Get a new molecule ready */
		if (atom->isHeteroAtom())
		{
			_myPolymer = PolymerPtr();
			_chainFrag = 0;
			std::string chain_id = atom->getChainID() + "0";

			/* Are we returning to a previous chain? */
			
			MoleculePtr testMol = _myCrystal->molecule(chain_id);
			
			if (testMol)
			{
				_myMolecule = testMol;
			}
			else
			{
				_myMolecule = MoleculePtr(new Molecule());
				_myMolecule->setAbsoluteBFacSubtract(0);
			}
		}
		else
		{
			_myPolymer = PolymerPtr(new Polymer());
			_myMolecule = _myPolymer;
		}

		if (_breakMe)
		{
			_chainFrag++;
			_breakMe = false;
		}
		else
		{
			_chainFrag = 0;
		}

		std::string chain_id = atom->getChainID() + i_to_str(_chainFrag);

		_myChain = atom->getChainID();
		_myMolecule->setChainID(chain_id);
		_myCrystal->addMolecule(_myMolecule);
	}
}

void PDBReader::validateResidue(AbsolutePtr atom)
{
	/* Add disordered placeholders! */

	if (!_myPolymer)
	{
		shout_at_helen("Something has gone wrong. The program\n"\
		               "is trying to ensure the residue is\n"\
		               "correct, but there is no protein to\n"\
		               "check!");
	}

	int difference = atom->getResNum() - _residueNum;

	/* Something crazy going on, warn user. */
	if (difference < 0 && _myMolecule->atomCount() != 0)
	{
		warn_user("Residue number negative or going backwards in PDB. "\
		          "Will try to cope. Will probably fail. "\
		"Residue " + i_to_str(atom->getResNum())
		+ " vs " + i_to_str(_residueNum));
	}
	else if (difference > 1 && _myMolecule->atomCount() != 0)
	{
		warn_user("Break in PDB chain " + atom->getChainID() + " of "\
		          + i_to_str(difference) + " residues.\n"\
		"Assigning to new PDB chain fragment. "\
		"Residue " + i_to_str(atom->getResNum()));

		_breakMe = true;
		validateMolecule(atom);
	}

	/* All set up already, no problem. */
	if (difference == 0 && _myMolecule->atomCount() > 0)
	{
		return;
	}

	/* Now we add the final residue as a new one */

	_residueNum = atom->getResNum();
	_myMonomer = MonomerPtr(new Monomer());
	_myMonomer->setResidueNum(_residueNum);
	_myMonomer->setIdentifier(atom->getResName());
	_myPolymer->addMonomer(_myMonomer);
	_myMonomer->setup();
}

void PDBReader::addAnisotropicBFactors(std::string line)
{
	if (line.length() < 70)
	{
		shout_at_user("This ANISOU line in the PDB\n"\
		              "is less than 70 characters and might be\n"\
		              "non-standard. Please check the formatting\n"\
		              "of this line.\n\n" + line);
	}

	std::string u11, u12, u13, u22, u23, u33;
	std::string resNum, chainID, atomName, resName;

	atomName = line.substr(11, 5);
	resName = line.substr(17, 3);
	chainID = line.substr(21, 2);
	resNum = line.substr(23, 7);
	u11 = line.substr(28, 7);
	u22 = line.substr(35, 7);
	u33 = line.substr(42, 7);
	u12 = line.substr(49, 7);
	u13 = line.substr(56, 7);
	u23 = line.substr(63, 7);

	double u11_val = atof(u11.c_str());
	double u12_val = atof(u12.c_str());
	double u13_val = atof(u13.c_str());
	double u22_val = atof(u22.c_str());
	double u23_val = atof(u23.c_str());
	double u33_val = atof(u33.c_str());

	const double scale = 1 / 10e3;
	mat3x3 tensor = make_mat3x3();
	tensor.vals[0] = u11_val * scale;
	tensor.vals[1] = u12_val * scale;
	tensor.vals[2] = u13_val * scale;
	tensor.vals[3] = u12_val * scale;
	tensor.vals[4] = u22_val * scale;
	tensor.vals[5] = u23_val * scale;
	tensor.vals[6] = u13_val * scale;
	tensor.vals[7] = u23_val * scale;
	tensor.vals[8] = u33_val * scale;

	if (_myAbsolute)
	{
		_myAbsolute->setTensor(tensor);
	}
	else
	{
		warn_user("Unassigned anisotropic B factors.");
	}
}

void PDBReader::addAtomToMolecule(std::string line)
{
	AbsolutePtr abs = makeAbsolute(line);
	/* Makes sure we're ready to append to the correct molecule */
	
	if (!abs)
	{
		return;
	}

	validateMolecule(abs);

	if (abs->isHeteroAtom())
	{
		if (abs->getResName() == "hoh")
		{
			abs->addToMolecule(_myHOH);
			abs->getAtom()->setWater(true);
		}
		else
		{
			abs->addToMolecule(_myMolecule);
		}
	}
	else
	{
		/* Ignore if hydrogen */
		if (abs->getElementSymbol() == "H")
		{
//			return;
		}

		/* Makes sure we're ready to append to the correct residue */
		validateResidue(abs);
		abs->addToMonomer(_myMonomer);
	}

}

void PDBReader::getSymmetry(std::string line)
{
	std::string aData, bData, cData, alphaData, betaData, gammaData, spaceGroup;

	aData = line.substr(6, 9);
	bData = line.substr(15, 9);
	cData = line.substr(24, 9);
	alphaData = line.substr(33, 7);
	betaData = line.substr(40, 7);
	gammaData = line.substr(47, 7);
	spaceGroup = line.substr(55, 10);

	double a = atof(aData.c_str());
	double b = atof(bData.c_str());
	double c = atof(cData.c_str());
	double alpha = atof(alphaData.c_str());
	double beta = atof(betaData.c_str());
	double gamma = atof(gammaData.c_str());

	CSym::CCP4SPG *spg = CSym::ccp4spg_load_by_ccp4_spgname(spaceGroup.c_str());
	_myCrystal->setUnitCell(a, b, c, alpha, beta, gamma);
	_myCrystal->setSpaceGroup(spg);
	_myCrystal->setupSymmetry();

	_foundCrystal = true;
}

void PDBReader::parseLine(std::string line)
{
	if (line.substr(0, 6) == "CRYST1")
	{
		getSymmetry(line);
		return;
	}

	if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")
	{
		addAtomToMolecule(line);
		return;
	}

	if (line.substr(0, 6) == "ANISOU")
	{
		addAnisotropicBFactors(line);
		return;
	}
	
	std::string rwork = "REMARK   3   R VALUE            (WORKING SET) : ";
	std::string rfree = "REMARK   3   FREE R VALUE                     : ";

	if (line.substr(0, rwork.size()) == rwork)
	{
		std::string rest = line.substr(rwork.size(), std::string::npos);
		_rWork = atof(rest.c_str());
	}
	else if (line.substr(0, rfree.size()) == rfree)
	{
		std::string rest = line.substr(rwork.size(), std::string::npos);
		_rFree = atof(rest.c_str());
	}
}

void PDBReader::parse()
{
	if (!file_exists(filename))
	{
		Shouter *shout;
		shout = new Shouter("PDB file " + filename + " does not exist.");
		throw shout;
	}
	
	/* Prepare water network for HOH atoms */
	_myHOH = WaterNetworkPtr(new WaterNetwork());
	_myHOH->setChainID("HOH");

	std::string pdbContents = get_file_contents(filename);

	std::vector<std::string> lines = split(pdbContents, '\n');

	for (int i = 0; i < lines.size(); i++)
	{
		parseLine(lines[i]);
	}
	
	if (_myHOH->atomCount())
	{
		_myCrystal->addMolecule(_myHOH);
	}
}

void PDBReader::setFilename(std::string file)
{
	filename = file;
}

CrystalPtr PDBReader::getCrystal()
{
	if (_myCrystal)
	{
		return _myCrystal;
	}

	std::cout << "Loading a crystal from PDB file " << filename << ".\n" << std::endl;

	_myCrystal = CrystalPtr(new Vagabond::Crystal());
	_myCrystal->setFilename(getBaseFilename(filename));
	Options::getRuntimeOptions()->addCrystal(_myCrystal);
	int count = Options::getRuntimeOptions()->crystalCount();
	Options::getRuntimeOptions()->setActiveCrystal(count - 1);

	_residueNum = -INT_MAX;
	_myChain = "";

	parse();

	if (!_foundCrystal)
	{
		shout_at_user("PDB file does not contain the CRYST1\n" \
		              "entry line which has details about the\n" \
		"unit cell dimensions and space group.");
	}

	return _myCrystal;
}

std::string PDBReader::writeCryst(mat3x3 frac2real, CSym::CCP4SPG *spg)
{
	double vals[6];
	unit_cell_from_mat3x3(frac2real, &vals[0]);

	std::ostringstream stream;
	stream << "CRYST1";
	stream << std::fixed;
	
	for (int i = 0; i < 3; i++)
	{
		stream << std::setw(9) << std::setprecision(3) << vals[i];
	}

	for (int i = 3; i < 6; i++)
	{
		stream << std::setw(7) << std::setprecision(2) << vals[i];
	}

	stream << " " << spg->symbol_xHM << std::endl;

	return stream.str();
}

std::string PDBReader::writeLine(AtomPtr atom, vec3 placement, int count,
                                 double occupancy, double bFactor)
{
	std::ostringstream stream;
	stream << atom->pdbLineBeginning();
	stream << std::fixed << std::setw(10) << std::setprecision(3) << placement.x;
	stream << std::fixed << std::setw(8) << std::setprecision(3) << placement.y;
	stream << std::fixed << std::setw(8) << std::setprecision(3) << placement.z;
	stream << std::fixed << std::setw(6) << std::setprecision(2) << occupancy;
	stream << std::fixed << std::setw(6) << std::setprecision(2) << bFactor;
	stream << "          ";
	stream << std::setw(2) << atom->getElement()->getSymbol();
	stream << "  " << std::endl;

	return stream.str();
}


