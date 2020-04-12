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
#include "FFT.h"
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
		_limit = -1;
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

	VagFFTPtr getFFT()
	{
		return _fft;
	}
	
	VagFFTPtr getOriginal()
	{
		return _original;
	}

	VagFFTPtr getPartial()
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
	
	void setResLimit(double limit)
	{
		_limit = limit;
	}

	void copyToFFT(VagFFTPtr vag);
protected:
	VagFFTPtr _fft;
	VagFFTPtr _partial;
	VagFFTPtr _original;
	std::string _filename;

	float _minRes;
	float _maxRes;
	float _limit;

private:
	/* Contains the reciprocal reflection values */

};

#endif /* defined(__vagabond__Diffraction__) */
