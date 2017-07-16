//
//  Model.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Model.h"
#include "shared_ptrs.h"
#include "Molecule.h"
#include "Atom.h"


void Model::addToMolecule(MoleculePtr molecule)
{
	molecule->addModel(shared_from_this());
}