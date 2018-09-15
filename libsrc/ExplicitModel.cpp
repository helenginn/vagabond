//
//  ExplicitModel.cpp
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "ExplicitModel.h"
#include "Options.h"
#include "fftw3d.h"
#include "CSV.h"
#include "Anisotropicator.h"

ExplicitModel::ExplicitModel()
{
	_modifySample = -1;
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


std::vector<vec3> ExplicitModel::polymerCorrectedPositions()
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
		
		if (positions->size() <= 1)
		{
			posOnly.push_back(subtract);
			break;
		}

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

std::vector<BondSample> ExplicitModel::getFinalPositions()
{
	if (!_recalcFinal && _finalSamples.size())
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
	
	
	/* Deset flag */
	
	_recalcFinal = false;
	
	if (_useMutex)
	{
		std::lock_guard<std::mutex> lock(guiLock);
		_finalPositions = posOnly;
		_absolute = meanOfManyPositions(&_finalSamples);
	}
	else
	{
		_finalPositions = posOnly;
		_absolute = meanOfManyPositions(&_finalSamples);
	}

	return _finalSamples;
}

void ExplicitModel::addRealSpacePositions(FFTPtr real, vec3 offset)
{
	std::vector<BondSample> positions = getFinalPositions();
	
	double realLimits = real->scales[0] * real->nx;
	vec3 absolute = getAbsolutePosition();

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 placement = positions[i].start;
		placement = vec3_add_vec3(placement, offset);
		vec3 relative = vec3_subtract_vec3(placement, absolute);

		double occupancy = positions[i].occupancy;

		vec3_mult(&relative, 1 / realLimits);

		if (relative.x != relative.x)
		{
			continue;
		}
		
		bool interpolate = false;
		if (interpolate)
		{
			FFT::collapseFrac(&relative.x, &relative.y, &relative.z);

			double x = relative.x * real->nx;
			double y = relative.y * real->ny;
			double z = relative.z * real->nz;
			
			real->addInterpolatedToReal(x, y, z, occupancy);
		}
		else
		{
			real->addBlurredToReal(relative.x, relative.y, relative.z,
			                       occupancy);
		}
	}
}

FFTPtr ExplicitModel::makeRealSpaceDistribution()
{
	double n = fftGridLength();
	vec3 offset = empty_vec3();
	
	if (_overrideN > 0) n = _overrideN;
	
	/* Don't panic, invert scale below... this is in real space */
	double maxDStar = Options::getRuntimeOptions()->getActiveCrystalDStar();
	double scale = 1.0 / (2 * maxDStar);

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);
	fft->createFFTWplan(1);

	addRealSpacePositions(fft, offset);

	fft->fft(1);
	fft->invertScale();

	FFTPtr newPtr;
	newPtr.reset(new FFT(*fft));
	return newPtr;
}

FFTPtr ExplicitModel::makeDistribution()
{
	return makeRealSpaceDistribution();
}

void ExplicitModel::writePositionsToFile(std::string filename)
{
	if (!hasExplicitPositions())
	{
		return;
	}

	CSVPtr csv = CSVPtr(new CSV(3, "x", "y", "z"));
	std::vector<BondSample> positions = getFinalPositions();
	
	for (int i = 0; i < positions.size(); i++)
	{
		csv->addEntry(3, positions[i].start.x, 
		              positions[i].start.y,
		              positions[i].start.z);
	}

	csv->writeToFile(filename);	
}

void ExplicitModel::refreshPositions()
{
	getFinalPositions();
}

void ExplicitModel::setPosN(int choice, double value)
{
	/*
	double *vec = &_position.x;

	if (_modifySample >= 0)
	{
		if (_modifySample >= _finalSamples.size())
		{
			return;
		}

		vec = &_bondSamples[_modifySample].start.x;
		_recalcFinal = true;
	}
	
	*(vec + choice) = value;
	*/
}

double ExplicitModel::getMeanSquareDeviation()
{
	getAnisotropy(true);
	return _isotropicAverage * 8 * M_PI * M_PI;
}


double ExplicitModel::getPosN(int choice)
{
	/*
	double *vec = &_position.x;
	if (_modifySample >= 0)
	{
		if (_modifySample >= _finalSamples.size())
		{
			return *(&_position.x + choice);
		}

		vec = &_bondSamples[_modifySample].start.x;
	}
	
	return *(vec + choice);
	*/
	return 0;
}

/* For GUI */
std::vector<vec3> ExplicitModel::fishPositions(vec3 *ave)
{
	std::vector<vec3> positions;
	std::lock_guard<std::mutex> lock(guiLock);
	positions = _finalPositions;
	if (ave)
	{
		*ave = _absolute;
	}

	return positions;
}

void ExplicitModel::getAnisotropy(bool withKabsch)
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
		_smallness = tropicator.smallness();
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
		_smallness = tropicator.smallness();
		_anisotropyExtent = tropicator.anisotropyExtent();
	}
}

vec3 ExplicitModel::longestAxis()
{
	getAnisotropy();
	return _longest;
}

double ExplicitModel::anisotropyExtent(bool withKabsch)
{
	getAnisotropy();
	return _anisotropyExtent;
}
