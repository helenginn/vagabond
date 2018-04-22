//
//  Bucket.h
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Bucket__
#define __vagabond__Bucket__

#include <stdio.h>
#include "shared_ptrs.h"
#include <vector>
#include "vec3.h"
#include "../libccp4/csymlib.h"
#include "../libccp4/ccp4_spg.h"

class Bucket
{
public:
	virtual void addSolvent() = 0;
	
	virtual ~Bucket() {}

	void scaleSolvent();
	static double scaleSolventScore(void *object);
	double scaleAndAddSolventScore();
	void applySymOps(CSym::CCP4SPG *spaceGroup, double res);
	void fourierTransform(int dir, double res);
	void writeMillersToFile(std::string prefix, double maxRes);
	void abandonCalculations();
	
	void setCrystal(CrystalPtr crystal)
	{
		_crystal = crystal;
	}
	
	void setData(DiffractionPtr data)
	{
		_data = data;
	}
	
	static double getSolvScale(void *object)
	{
		return static_cast<Bucket *>(object)->_solvScale;
	}
	
	static void setSolvScale(void *object, double value)
	{
		static_cast<Bucket *>(object)->_solvScale = value;
	}

	static double getSolvBFac(void *object)
	{
		return static_cast<Bucket *>(object)->_solvBFac;
	}
	
	static void setSolvBFac(void *object, double value)
	{
		static_cast<Bucket *>(object)->_solvBFac = value;
	}

protected:
	CrystalWkr _crystal;
	FFTPtr _solvent;
	DiffractionPtr _data;
	
	CrystalPtr getCrystal()
	{
		return _crystal.lock();
	}
private:
	double _solvScale;
	double _solvBFac;
};

#endif /* defined(__vagabond__Bucket__) */
