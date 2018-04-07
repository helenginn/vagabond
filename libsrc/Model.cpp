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
#include "Anisotropicator.h"
#include "Element.h"

bool Model::_useMutex = false;

Model::Model()
{
	_recalcFinal = true;
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
	double scale = 2 * MAX_SCATTERING_DSTAR;
	/* Target dimension in Angstroms */
	double dim = biggestStdevDim() * 2;
	
	int some = 4;
	
	if (getAtom()->getElement()->electronCount() <= 1)
	{
		some = 2;	
	}
	
	/* Add some Angstroms for good luck */
	dim += some;
	dim *= 2;
	
	int n = scale * dim + 0.5;
	//n = ATOM_SAMPLING_COUNT;
	return n;
}

std::vector<vec3> Model::polymerCorrectedPositions()
{
	std::vector<vec3> posOnly;
	std::vector<BondSample> *positions = getManyPositions();
	posOnly.reserve(positions->size());

	MoleculePtr molecule = getMolecule();
	std::vector<vec3> offsets;
	std::vector<vec3> transTensorOffsets;
	std::vector<vec3> rotationCentres;
	std::vector<mat3x3> rotations;
	std::vector<mat3x3> extraRotations;
	vec3 rotCentre;

	if (molecule)
	{
		offsets = molecule->getCentroidOffsets();
		rotations = molecule->getRotationCorrections();
		rotationCentres = molecule->getRotationCentres();
		transTensorOffsets = molecule->getTransTensorOffsets();
		extraRotations = molecule->getExtraRotations();
		rotCentre = molecule->getExtraRotationCentre();
	}

	bool hasCentre = true;

	if (rotCentre.x != rotCentre.x)
	{
		hasCentre = false;
	}

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 subtract = positions->at(i).start;

		// perform rotation element of superposition

		if (rotations.size() > i && rotationCentres.size() > i)
		{
			vec3 tmp = vec3_add_vec3(subtract, rotationCentres[i]);
			mat3x3_mult_vec(rotations[i], &tmp);
			subtract = vec3_subtract_vec3(tmp, rotationCentres[i]);
		}

		// remove translation aspect of superposition

		if (offsets.size() > i)
		{
			subtract = vec3_subtract_vec3(subtract, offsets[i]);
		}

		// additional rotations applying to whole-molecule

		if (extraRotations.size() > i && rotationCentres.size() > i)
		{
			vec3 *bestVec = &rotationCentres[i];

			if (hasCentre)
			{
				bestVec = &rotCentre;
			}

			vec3 tmp = vec3_add_vec3(subtract, *bestVec); 
			mat3x3_mult_vec(extraRotations[i], &tmp);
			subtract = vec3_subtract_vec3(tmp, *bestVec);
		}

		// remove the translation from the tensor for translation

		if (transTensorOffsets.size() > i)
		{
			subtract = vec3_subtract_vec3(subtract, transTensorOffsets[i]);
		}

		posOnly.push_back(subtract);
	}

	return posOnly;
}

vec3 meanOfManyPositions(std::vector<BondSample> *positions)
{
	vec3 sum = make_vec3(0, 0, 0);
	double weights = 0;

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 tmp = (*positions)[i].start;
		vec3_mult(&tmp, (*positions)[i].occupancy);
		sum = vec3_add_vec3(sum, tmp);
		weights += (*positions)[i].occupancy;
	}

	vec3_mult(&sum, 1 / weights);

	return sum;
}

std::vector<BondSample> Model::getFinalPositions()
{
	if (!_recalcFinal)
	{
		return _finalSamples;	
	}
	
	/* Recalculate final positions */

	std::vector<BondSample> *positions = getManyPositions();
	_finalSamples = *positions;
	std::vector<vec3> posOnly = polymerCorrectedPositions();

	for (int i = 0; i < _finalSamples.size(); i++)
	{
		if (i < posOnly.size())
		{
			_finalSamples[i].start = posOnly[i];
		}
	}

	_absolute = meanOfManyPositions(&_finalSamples);
	
	/* Deset flag */
	
	_recalcFinal = false;
	
	if (_useMutex)
	{
		std::lock_guard<std::mutex> lock(guiLock);
		_finalPositions = posOnly;
	}
	else
	{
		_finalPositions = posOnly;
	}

	return _finalSamples;
}

void Model::propagateChange(int depth, bool refresh)
{
	_recalcFinal = true;
}


/* For GUI */
std::vector<vec3> Model::fishPositions()
{
	std::vector<vec3> positions;
	std::lock_guard<std::mutex> lock(guiLock);
	positions = _finalPositions;

	return positions;
}


void Model::getAnisotropy(bool withKabsch)
{
	if (withKabsch)
	{
		std::vector<BondSample> finals = getFinalPositions();
		std::vector<vec3> finalPoints;

		for (int i = 0; i < finals.size(); i++)
		{
			finalPoints.push_back(finals[i].start);
		}

		if (!finalPoints.size()) return;
		
		Anisotropicator tropicator;
		tropicator.setPoints(finalPoints);
		_realSpaceTensor = tropicator.getTensor();
		_anisotropyExtent = tropicator.anisotropyExtent();
		_longest = tropicator.longestAxis();
		_isotropicAverage = tropicator.isotropicAverage();
	}
	else
	{
		std::vector<BondSample> *positions = getManyPositions();
		std::vector<vec3> points;

		for (int i = 0; i < positions->size(); i++)
		{
			points.push_back((*positions)[i].start);
		}

		Anisotropicator tropicator;
		tropicator.setPoints(points);
		_anisotropyExtent = tropicator.anisotropyExtent();
	}
}

vec3 Model::longestAxis()
{
	getAnisotropy(true);
	return _longest;
}

double Model::anisotropyExtent(bool withKabsch)
{
	getAnisotropy(withKabsch);
	return _anisotropyExtent;
}


mat3x3 Model::getRealSpaceTensor()
{
	longestAxis();
	return _realSpaceTensor;
}

void Model::addProperties()
{
	addMat3x3Property("real_space_tensor", &_realSpaceTensor);
}
