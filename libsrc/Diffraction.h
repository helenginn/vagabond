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

class Diffraction : public Dataset
{
public:
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
protected:
	FFTPtr fft;
	std::string _filename;

private:
	/* Contains the reciprocal reflection values */

};

#endif /* defined(__vagabond__Diffraction__) */
