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

class PDBReader
{
public:

	PDBReader();
	void setFilename(std::string);
	CrystalPtr getCrystal();

private:
	std::string filename;

	/* Keep track of what we're dealing with */
	MoleculePtr _myMolecule;
	PolymerPtr _myPolymer;
	CrystalPtr _myCrystal;
	MonomerPtr _myMonomer;
	std::string _myChain;
	int _residueNum;

	void getSymmetry(std::string line);
	AbsolutePtr makeAtom(std::string line);
	void parseLine(std::string line);
	void parse();
	void validateMolecule(AbsolutePtr atom);
	void validateResidue(AbsolutePtr atom);
	void addAtomToMolecule(std::string line);

	bool _foundCrystal;
};

#endif /* defined(__vagabond__PDBReader__) */
