//
//  Model.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Model.h"
#include "CSV.h"
#include "shared_ptrs.h"
#include "Polymer.h"
#include "Monomer.h"
#include "Atom.h"
#include "Element.h"
#include "Shouter.h"

bool Model::_useMutex = false;

Model::Model()
{
	_recalcFinal = true;
	_recalcDist = true;
	_realSpaceTensor = make_mat3x3();
}

void Model::addToMonomer(MonomerPtr monomer)
{
	monomer->addModel(shared_from_this());

	PolymerPtr polymer = monomer->getPolymer();
	MoleculePtr molecule = boost::static_pointer_cast<Molecule>(polymer);
}

double returnOne(void *object, double x, double y, double z)
{
	return 1;
}

double Model::biggestStdevDim()
{
	getAnisotropy(true);
	double max = 0;
	mat3x3 mat = getRealSpaceTensor();
	
	for (int i = 0; i < 3; i++)
	{
		vec3 axis = mat3x3_axis(mat, i);
		double sqlength = vec3_sqlength(axis);
		
		if (sqlength > max) max = sqlength;
	}
	
	return sqrt(max);
}

int Model::fftGridLength()
{
	double maxDStar = Options::getRuntimeOptions()->getActiveCrystalDStar();
	double scale = 2.0 * maxDStar;

	/* Target dimension in Angstroms */
	double dim = biggestStdevDim() * 2;
	
	if (dim != dim || dim <= 0)
	{
		dim = 0;
	}
	
	/* Add some Angstroms for good luck */
	dim += 2;
	dim *= 2;
	
	int n = scale * dim + 0.5;
	
	if (n % 2 == 1)
	{
		n += 1;
	}
	
	return n;
}

FFTPtr Model::getDistribution()
{
	if (_recalcDist)
	{
		_lastDistribution = makeDistribution();
	}
	
	return _lastDistribution;
}

void Model::propagateChange(int depth, bool refresh)
{
	_recalcFinal = true;
}

vec3 Model::longestAxis()
{
	return _longest;
}

double Model::anisotropyExtent(bool withKabsch)
{
	return _anisotropyExtent;
}

double Model::smallness()
{
	return _smallness;
}

mat3x3 Model::getRealSpaceTensor()
{
	longestAxis();
	return _realSpaceTensor;
}

void Model::addProperties()
{
	addMat3x3Property("real_space_tensor", &_realSpaceTensor, true);
}




