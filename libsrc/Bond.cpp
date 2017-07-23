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

Bond::Bond(AtomPtr major, AtomPtr minor)
{
	_major = major;
	_minor = minor;
}

void Bond::setup()
{
	getMinor()->setModel(shared_from_this());
}

double Bond::getVoxelValue(void *obj, double x, double y, double z)
{
	return 1;
}

FFTPtr Bond::getDistribution()
{
	return FFTPtr();// prepareDistribution(&getVoxelValue, 1);
}