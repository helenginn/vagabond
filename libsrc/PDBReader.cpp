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
#include "Molecule.h"
#include "Absolute.h"
#include "shared_ptrs.h"
#include <vector>
#include "Shouter.h"
#include "mat3x3.h"
#include "Crystal.h"

PDBReader::PDBReader()
{
	_foundCrystal = false;
}

void PDBReader::addAtom(std::string line)
{
	std::string xData, yData, zData, element, bFactor, occupancy;

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

	vec3 vec = make_vec3(xValue, yValue, zValue);
	AbsolutePtr abs = AbsolutePtr(new Absolute(vec, bFacValue, element, occValue));
	abs->addToMolecule(myMolecule);
}

void PDBReader::getSymmetry(std::string line)
{
	std::string aData, bData, cData, alphaData, betaData, gammaData, spaceGroup;

	if (line.length() < 65)
	{
		printf("Error.\n");
		exit(1);
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

	myCrystal->setHKL2Real(hkl2real);
	myCrystal->setReal2HKL(real2hkl);

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
		addAtom(line);
	}
}

void PDBReader::parse()
{
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
	myCrystal = CrystalPtr(new Crystal());
	myMolecule = MoleculePtr(new Molecule());

	parse();

	if (!_foundCrystal)
	{
		shout_at_user("PDB file does not contain the CRYST1\n" \
					  "entry line which has details about the\n" \
					  "unit cell dimensions and space group.");
	}

	myCrystal->addMolecule(myMolecule);

	std::cout << "Atoms: " << myMolecule->atomCount() << std::endl;

	return myCrystal;
}
