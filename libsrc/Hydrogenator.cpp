//
// Hydrogenator.cpp
// vagabond
//
// Created by Helen Ginn on 22/03/2018
// Copyright (c) 2018 Helen Ginn
//

#include "Hydrogenator.h"

Hydrogenator::Hydrogenator()
{
	
}


void Hydrogenator::hydrogenate()
{
	if (!_monomer)
	{
		return;
	}
	
	BackbonePtr bone = _monomer->getBackbone();

	if (!bone)
	{
		return;	
	}
	
	AtomPtr nitrogen = bone->findAtom("N");
	ElementPtr hydrogenElement = Elemen 
	
	if (nitrogen)
	{
		AtomPtr nHydrogen = AtomPtr(new Atom());
		nHydrogen->setFromPDB(false);

		nHydrogen->setInitialBFactor(nitrogen->getInitialBFactor());
		nHydrogen->setElement(hydrogenElement);
		nHydrogen->setAtomName("H");
		nHydrogen->setOriginalOccupancy(1.);
	}
}
