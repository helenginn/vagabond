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
#include "Anisotropicator.h"
#include "Element.h"
#include "Shouter.h"
#include "fftw3d.h"

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
	
	int some = 2;

	if (getAtom()->getElement()->electronCount() <= 1)
	{
		some = 2;	
	}
	
	/* Add some Angstroms for good luck */
	dim += some;
	dim *= 2;
	
	int n = scale * dim + 0.5;
	
	if (n % 2 == 1)
	{
		n += 1;
	}
	
//	std::cout << "Size " << dim << ", scale " << scale << ", n = " << n << std::endl;

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

double Model::smallness()
{
	getAnisotropy(true);
	return _smallness;
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


void Model::addRealSpacePositions(FFTPtr real, vec3 offset)
{
	if (!hasExplicitPositions())
	{
		shout_at_helen("Inexplicit model asked to add positions");	
	}
	
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

FFTPtr Model::makeRealSpaceDistribution()
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

void Model::writePositionsToFile(std::string filename)
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

