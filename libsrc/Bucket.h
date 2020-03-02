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
#include "PartialStructure.h"
#include <vector>
#include <map>
#include "vec3.h"
#include "../libccp4/csymlib.h"
#include "../libccp4/ccp4_spg.h"


#define MAX_CHECK_DISTANCE 6.0
#define MIN_CHECK_DISTANCE 1.0


/**
 * \class Bucket
 * \brief Abstract class for adding solvent (via implemented method),
 * applying symmetry operations and scaling.
 */

class Bucket : public PartialStructure
{
public:
	Bucket();
	virtual void addSolvent() = 0;
	virtual void reportScale();
	
	static BucketPtr chosenBucket();
	
	virtual ~Bucket()
	{
		_atomPtrs.clear();
		std::vector<Atom *>().swap(_atomPtrs);
	}

	int getRandomValues(double left, double *right, double *angle);
	int getReallyRandomValues(double left, double *right, double *angle);

	void fourierTransform(int dir);
	void writeMillersToFile(std::string prefix, double maxRes);
	void abandonCalculations();
	
	virtual void postScaleWork() {};
	
	/* only use before FFT */

	bool isSolvent(int index);
	Atom *nearbyAtom(int index);
	
	VagFFTPtr getSolvent()
	{
		return _solvent;
	}
protected:
	VagFFTPtr _solvent;
	std::vector<Atom *> _atomPtrs;
	
private:
	/* Mask regions with protein = 0, solvent = 1 and protein/solvent
	* interface = 2 */
	
	double _averages[3];
};

#endif /* defined(__vagabond__Bucket__) */
