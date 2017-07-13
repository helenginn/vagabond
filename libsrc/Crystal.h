//
//  Crystal.h
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Crystal__
#define __vagabond__Crystal__

#include <stdio.h>
#include <vector>
#include "shared_ptrs.h"
#include "mat3x3.h"

class Crystal
{
public:
	void addMolecule(MoleculePtr molecule)
	{
		_molecules.push_back(molecule);
	}

	void setReal2HKL(mat3x3 mat);
	void setHKL2Real(mat3x3 mat);


private:
	std::vector<MoleculePtr> _molecules;

	mat3x3 _hkl2real;
	mat3x3 _real2hkl;

};

#endif /* defined(__vagabond__Crystal__) */
