//
//  PDBReader.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__PDBReader__
#define __vagabond__PDBReader__

#include <stdio.h>
#include <iostream>
#include "shared_ptrs.h"
#include "vec3.h"

class PDBReader
{
public:

	PDBReader();
	void setFilename(std::string);
	CrystalPtr getCrystal();

	static std::string writeLine(AtomPtr atom, vec3 placement, int count,
								 double occupancy, double bFactor);

private:
	std::string filename;

	/* Keep track of what we're dealing with */
	MoleculePtr _myMolecule;
	PolymerPtr _myPolymer;
	CrystalPtr _myCrystal;
	MonomerPtr _myMonomer;
	std::string _myChain;
	bool _breakMe;
	int _chainFrag;
	AbsolutePtr _myAbsolute;
	int _residueNum;

	std::string longChainID(std::string base);

	void getSymmetry(std::string line);
	AbsolutePtr makeAbsolute(std::string line);
	void parseLine(std::string line);
	void parse();
	void validateMolecule(AbsolutePtr atom);
	void validateResidue(AbsolutePtr atom);
	void addAtomToMolecule(std::string line);
	void addAnisotropicBFactors(std::string line);

	bool _foundCrystal;
};

#endif /* defined(__vagabond__PDBReader__) */
