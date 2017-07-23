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

PDBReader::PDBReader()
{
	_foundCrystal = false;
}

AbsolutePtr PDBReader::makeAtom(std::string line)
{
	if (line.length() < 78)
	{
		shout_at_user("This ATOM line in the PDB\n"\
					  "is less than 80 characters and must be\n"\
					  "non-standard. Please check the formatting\n"\
					  "of this line.\n\n" + line);
	}

	std::string xData, yData, zData, element, bFactor, occupancy;
	std::string resNum, chainID, atomName, resName;

	atomName = line.substr(13, 4);
	resName = line.substr(17, 3);
	chainID = line.substr(21, 2);
	resNum = line.substr(23, 7);
	xData = line.substr(30, 8);
	yData = line.substr(38, 8);
	zData = line.substr(46, 8);
	occupancy = line.substr(54, 6);
	bFactor = line.substr(60, 6);
	element = line.substr(77, 1);

	double xValue = atof(xData.c_str());
	double yValue = atof(yData.c_str());
	double zValue = atof(zData.c_str());
	double bFacValue = atof(bFactor.c_str());
	double occValue = atof(occupancy.c_str());
	double resNumValue = atoi(resNum.c_str());

	vec3 vec = make_vec3(xValue, yValue, zValue);
	AbsolutePtr abs = AbsolutePtr(new Absolute(vec, bFacValue, element, occValue));
	abs->setIdentity(resNumValue, chainID, resName, atomName);

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

			if (_myMolecule->className() == "Polymer")
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
				  "Will try to cope. Will probably fail.");
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
	_myMonomer->setup();
	_myMonomer->setResidueNum(_residueNum);
	_myMonomer->setIdentifier(atom->getResName());
	_myPolymer->addMonomer(_myMonomer);

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
}

void PDBReader::parse()
{
	if (!FileReader::exists(filename))
	{
		shout_at_user("Cannot open file " + filename + ".");
	}

	std::string pdbContents = FileReader::get_file_contents(filename.c_str());

	std::vector<std::string> lines = FileReader::split(pdbContents, '\n');

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

	_myCrystal->tieAtomsUp();

	return _myCrystal;
}
