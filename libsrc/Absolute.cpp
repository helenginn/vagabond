//
//  Absolute.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "shared_ptrs.h"
#include "Absolute.h"
#include "Molecule.h"
#include "Atom.h"

Absolute::Absolute(vec3 pos, double bFac)
{
	position = pos;
	bFactor = bFac;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
	molecule->addModel(shared_from_this());

	AtomPtr myAtom = AtomPtr(new Atom());
	myAtom->addConnection(shared_from_this());

	molecule->addAtom(myAtom);
	myAtom->setPosition(position);
	_atom = myAtom;
}