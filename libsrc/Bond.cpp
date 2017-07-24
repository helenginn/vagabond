//
//  Bond.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bond.h"
#include "Atom.h"
#include "fftw3d.h"
#include <iostream>

Bond::Bond(AtomPtr major, AtomPtr minor)
{
	_major = major;
	_minor = minor;
	_torsionBasis = make_mat3x3();
	_torsionRadians = 0;
	_bondLength = 0;
	_minorGeometry = BondGeometryTetrahedral;

	vec3 majorPos = getMajor()->getPosition();
	vec3 minorPos = getMinor()->getPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	_bondLength = vec3_length(difference);
}

void Bond::setup()
{
	getMinor()->setModel(shared_from_this());
}

/* Do I really need this? */
void Bond::addDownstreamAtom(AtomPtr down)
{
	_downstreamAtoms.push_back(down);
}

void Bond::setAlignmentAtoms(AtomPtr heavyAlign, AtomPtr lightAlign)
{
	_heavyAlign = heavyAlign;
	_lightAlign = lightAlign;
	
}

double Bond::getVoxelValue(void *obj, double x, double y, double z)
{
	double distSq = x * x + y * y + z * z;

	if (sqrt(distSq) < 0.16)
	{
		return 1;
	}

	return 0;
}

FFTPtr Bond::getDistribution()
{
	/* Getting the distribution for the 'major' atom */

	FFTPtr inherited = getMajor()->getBlur();

	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2.0 * MAX_SCATTERING_DSTAR);

	/* Deal with a needed width of 4.0 Å for now */

	prepareDistribution(n, scale, this, &getVoxelValue);
	FFTPtr mine = getDistributionCopy();
	mine->fft(1);
	FFT::multiply(mine, inherited);

	return mine;
}