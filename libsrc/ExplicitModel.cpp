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
