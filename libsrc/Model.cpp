//
//  Model.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Model.h"
#include "shared_ptrs.h"
#include "Polymer.h"
#include "Monomer.h"
#include "Atom.h"

Model::Model()
{
	_realSpaceTensor = make_mat3x3();
}

void Model::addToMolecule(MoleculePtr molecule)
{
	molecule->addModel(shared_from_this());
}

void Model::addToMonomer(MonomerPtr monomer)
{
	monomer->addModel(shared_from_this());

	PolymerPtr polymer = monomer->getPolymer();
	MoleculePtr molecule = std::static_pointer_cast<Molecule>(polymer);
	molecule->addModel(shared_from_this());
}

double returnOne(void *object, double x, double y, double z)
{
	return 1;
}

FFTPtr Model::getZeroDistribution()
{
	double n = ATOM_SAMPLING_COUNT;
	double scale = 2 * MAX_SCATTERING_DSTAR;
	prepareDistribution(n, scale, this, &returnOne);

	return getDistributionCopy();
}

std::vector<BondSample> Model::getFinalPositions()
{
	std::vector<BondSample> *positions = getManyPositions(BondSampleThorough);
	std::vector<BondSample> copyPos = *positions;

	return copyPos;
}

void Model::propagateChange(bool activeGroupOnly)
{
	
}
