//
//  Diffraction.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Diffraction__
#define __vagabond__Diffraction__

#include <stdio.h>
#include "Dataset.h"
#include "fftw3d.h"
#include "shared_ptrs.h"
#include <string>

/**
 * \class Diffraction
 * \brief Dataset which provides diffraction intensities of some description.
 */

class Diffraction : public Dataset
{
public:
	Diffraction()
	{
		_minRes = 0;
		_maxRes = 0;
	}

	virtual ~Diffraction() {};

	void setFilename(std::string filename)
	{
		_filename = filename;
	}

	std::string getFilename()
	{
		return _filename;
	}

	virtual void load() = 0;

	FFTPtr getFFT()
	{
		return fft;
	}

	FFTPtr getPartial()
	{
		return _partial;
	}
	
	bool hasPartial()
	{
		return (!_partial);
	}

	double getMaxResolution()
	{
		return _maxRes;
	}
protected:
	FFTPtr fft;
	FFTPtr _partial;
	std::string _filename;

	float _minRes;
	float _maxRes;

private:
	/* Contains the reciprocal reflection values */

};

#endif /* defined(__vagabond__Diffraction__) */
