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
	MoleculePtr getMolecule();

private:
	std::string filename;
	MoleculePtr myMolecule;

	void getSymmetry(std::string line);
	void addAtom(std::string line);
	void parseLine(std::string line);
	void parse();

	bool _foundCrystal;
};

#endif /* defined(__vagabond__PDBReader__) */
