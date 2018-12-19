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

vec3 ExplicitModel::meanOfManyPositions(std::vector<BondSample> *positions)
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
	
	if (weights <= 0)
	{
		weights = 1;
	}

	vec3_mult(&sum, 1 / weights);

	return sum;
}

void ExplicitModel::sanityCheck()
{
	for (int i = 0; i < _storedSamples.size(); i++)
	{
		vec3 start = _storedSamples[i].start;
		vec3 old_start = _storedSamples[i].old_start;
		mat3x3 basis = _storedSamples[i].basis;
		
		if (start.x != start.x) 
		{
			std::cout << "Start position is nan for " << shortDesc() <<
			std::endl;
			return;
		}
		if (old_start.x != old_start.x)
		{
			std::cout << "Old start position is nan for " << shortDesc() <<
			std::endl;
			return;
		}
		
		for (int j = 0; j < 9; j++)
		{
			if (basis.vals[j] != basis.vals[j])
			{
				std::cout << "Basis value " << j << " is nan for " <<
				shortDesc() << std::endl;
				return;
			}
		}
	}
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
	
	std::vector<vec3> posOnly;
	posOnly.reserve(_finalSamples.size());

	for (int i = 0; i < _finalSamples.size(); i++)
	{
		posOnly.push_back(_finalSamples[i].start);
	}
	
	/* Deset flag to allow caching */
	
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

mat3x3 ExplicitModel::makeTorsionBasis(vec3 hPos, vec3 maPos,
                              vec3 miPos, vec3 lPos, double *newAngle)
{
	vec3 a2p = vec3_subtract_vec3(hPos, maPos);
	vec3 a2l = vec3_subtract_vec3(lPos, maPos);
	vec3 a2b = vec3_subtract_vec3(miPos, maPos);

	double sql = vec3_sqlength(a2b);
	double dotp = vec3_dot_vec3(a2p, a2b);
	double dotl = vec3_dot_vec3(a2l, a2b);
	double distp = dotp / sql;
	double distl = dotl / sql;

	vec3 a2bp = a2b;
	vec3 a2bl = a2b;
	vec3_mult(&a2bp, distp);
	vec3_mult(&a2bl, distl);
	vec3 heavy_join = vec3_add_vec3(maPos, a2bp);
	vec3 light_join = vec3_add_vec3(maPos, a2bl);

	vec3 reverse_bond = vec3_subtract_vec3(miPos, maPos);
	vec3 xNew = vec3_subtract_vec3(hPos, heavy_join);
	mat3x3 basis = mat3x3_rhbasis(xNew, reverse_bond);

	if (newAngle)
	{
		vec3 light = vec3_subtract_vec3(lPos, light_join);
		double angle = vec3_angle_with_vec3(light, xNew);

		vec3 recross = vec3_cross_vec3(light, xNew);
		double cosine = vec3_cosine_with_vec3(recross, reverse_bond);

		if (cosine > 0)
		{
			angle *= -1;
		}

		*newAngle = angle;
	}

	return basis;
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

	CSVPtr csv = CSVPtr(new CSV(4, "x", "y", "z", "k"));
	std::vector<BondSample> positions = getFinalPositions();
	
	for (int i = 0; i < positions.size(); i++)
	{
		csv->addEntry(3, positions[i].start.x, 
		              positions[i].start.y,
		              positions[i].start.z,
		              positions[i].kickValue);
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
