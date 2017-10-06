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
#include "Monomer.h"
#include "Absolute.h"
#include "csymlib.h"

PDBReader::PDBReader()
{
	_foundCrystal = false;
}

AbsolutePtr PDBReader::makeAtom(std::string line)
{
	if (line.length() < 78)
	{
		shout_at_user("This ATOM line in the PDB\n"\
					  "is less than 78 characters and might be\n"\
					  "non-standard. Please check the formatting\n"\
					  "of this line.\n\n" + line);
	}

	std::string xData, yData, zData, element, bFactor, occupancy;
	std::string resNum, chainID, atomName, atomNum, resName;

	atomNum = line.substr(6, 5);
	atomName = line.substr(11, 5);
	resName = line.substr(17, 3);
	chainID = line.substr(21, 2);
	resNum = line.substr(23, 7);
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
	double resNumValue = atoi(resNum.c_str());
	double atomNumValue = atoi(atomNum.c_str());

	vec3 vec = make_vec3(xValue, yValue, zValue);
	AbsolutePtr abs = AbsolutePtr(new Absolute(vec, bFacValue, element, occValue));
	_myAbsolute = abs;
	abs->setIdentity(resNumValue, chainID, resName, atomName, atomNumValue);

	if (line.substr(0, 6) == "HETATM")
	{
		abs->setHeteroAtom(true);
	}

	return abs;
}

void PDBReader::validateMolecule(AbsolutePtr atom)
{
	std::string newChain = atom->getChainID();

	if (newChain != _myChain)
	{
		/* We have switched proteins. Have we gone back to an old one? */

		if (_myCrystal->molecule(newChain))
		{
			_myMolecule = _myCrystal->molecule(newChain);

			if (_myMolecule->getClassName() == "Polymer")
			{
				_myPolymer = std::static_pointer_cast<Polymer>(_myMolecule);
			}
		}
		else
		{
			if (atom->isHeteroAtom())
			{
				_myPolymer = PolymerPtr();
				_myMolecule = MoleculePtr(new Molecule());
			}
			else
			{
				_myPolymer = PolymerPtr(new Polymer());
				_myMolecule = _myPolymer;
			}

			_myMolecule->setChainID(newChain);
			_myCrystal->addMolecule(_myMolecule);
		}
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
	if (difference < 0)
	{
		warn_user("Residue number going backwards in PDB. "\
				  "Will try to cope. Will probably fail."\
				  "residue " + i_to_str(atom->getResNum())
				  + " vs " + i_to_str(_residueNum));
	}

	/* All set up already, no problem. */
	if (difference == 0)
	{
		return;
	}

	/* More than one residue difference? Must be missing some
	 * disordered residues. We care. Add placeholders. */
	if (difference > 1)
	{
		difference--;
		_myPolymer->addUnknownMonomers(difference);
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
	AbsolutePtr abs = makeAtom(line);
	validateMolecule(abs);

	if (abs->isHeteroAtom())
	{
		abs->addToMolecule(_myMolecule);
	}
	else
	{
		validateResidue(abs);
		abs->addToMonomer(_myMonomer);
	}

}

void PDBReader::getSymmetry(std::string line)
{
	std::string aData, bData, cData, alphaData, betaData, gammaData, spaceGroup;

	if (line.length() < 65)
	{
		shout_at_user("The symmetry line (CRYST1) in the PDB\n"\
					  "is less than 65 characters and must be\n"\
					  "non-standard. Please check the formatting\n"\
					  "of this line.");
	}

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

	// finish me
	mat3x3 hkl2real = mat3x3_from_unit_cell(a, b, c, alpha, beta, gamma);
	mat3x3 real2hkl = mat3x3_inverse(hkl2real);

	CSym::CCP4SPG *spg = CSym::ccp4spg_load_by_spgname(spaceGroup.c_str());

	_myCrystal->setUnitCell(a, b, c, alpha, beta, gamma);
	_myCrystal->setSpaceGroup(spg);
	_myCrystal->setHKL2Real(hkl2real);
	_myCrystal->setReal2HKL(real2hkl);

	_foundCrystal = true;
}

void PDBReader::parseLine(std::string line)
{
	if (line.substr(0, 6) == "CRYST1")
	{
		getSymmetry(line);
	}

	if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")
	{
		addAtomToMolecule(line);
	}

	if (line.substr(0, 6) == "ANISOU")
	{
		addAnisotropicBFactors(line);
	}
}

void PDBReader::parse()
{
	if (!file_exists(filename))
	{
		shout_at_user("Cannot open file " + filename + ".");
	}

	std::string pdbContents = get_file_contents(filename);

	std::vector<std::string> lines = split(pdbContents, '\n');

	for (int i = 0; i < lines.size(); i++)
	{
		parseLine(lines[i]);
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

	_myCrystal = CrystalPtr(new Crystal());
	_myCrystal->setFilename(getBaseFilename(filename));

	_residueNum = 0;
	_myChain = "";

	parse();

	_myCrystal->summary();
	
	if (!_foundCrystal)
	{
		shout_at_user("PDB file does not contain the CRYST1\n" \
					  "entry line which has details about the\n" \
					  "unit cell dimensions and space group.");
	}

	return _myCrystal;
}
