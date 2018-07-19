//
//  Distributor.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Distributor__
#define __vagabond__Distributor__

#include <stdio.h>
#include "shared_ptrs.h"
#include "fftw3d.h"
#include <map>

typedef double(get_voxel_value)(void *obj, double x, double y, double z);

/**
* \class Distributor
* \brief Creates FFT grids for reciprocal space distributions.
*
* A number of classes require population of a three-dimensional grid
* containing reciprocal space information in preparation for a Fourier
* transform. This is an abstracted way to make sure that these classes provide
* suitable transforms, and re-use old calculations where necessary.
*
* A notable exception to the subclassing of Distributor is Bond. This is
* because Bond-derived distributions are calculated in real space.
*
*/

class Distributor
{
public:
	Distributor()
	{
		_calculated = false;
		_activeNum = 0;
		_overrideN = -1;
	}

	/** Call this if the subclass has changed (e.g., an Absolute atom has had
	* its B factor altered), and an appropriate FFT will be recalculated on
	* demand. */
	void recalculate()
	{
		_precalcFFTs.clear();
	}

	void overrideLength(int n)
	{
		_overrideN = n;
	}
protected:
	bool _calculated;
	int _overrideN;

	/** Get a copy of the distribution which is then alterable without
	* affecting the stored version. */
	FFTPtr getDistributionCopy()
	{
		if (_activeNum <= 0) return FFTPtr();

		FFTPtr newPtr;
		newPtr.reset(new FFT(*_precalcFFTs[_activeNum]));
		return newPtr;
	}

	/** Prepare a distribution if it has not already been made for a given
	* size.
	* \param n number of voxels along each edge
	* \param scale reciprocal scale in Ang^-1.
	* \param object pointer to an object to be called upon when filling the FFT with values.
	* \param voxel_value function to be called in order to fill the FFT with values, also passing in the object parameter.
	* \return reciprocal space FFT.
	*/
	FFTPtr prepareDistribution(double n, double scale, void *object,
	                           get_voxel_value *voxel_value);

	FFTPtr getFFT(int n)
	{
		return _precalcFFTs[n];	
	}
	
	virtual std::string getClassName()
	{
		return "Distributor";
	}

private:
	int _activeNum;
	std::map<int, FFTPtr> _precalcFFTs;
};

#endif /* defined(__vagabond__Distributor__) */
