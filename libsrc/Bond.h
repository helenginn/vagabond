//
//  Bond.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Bond__
#define __vagabond__Bond__

#include <stdio.h>
#include "Model.h"
#include <vector>
#include "mat3x3.h"
#include "Distributor.h"

typedef enum
{
	BondGeometryNone,
	BondGeometryTetrahedral,
} BondGeometryType;

class Bond : public Model, public Distributor
{
public:
	Bond(AtomPtr major, AtomPtr minor);
	void setup();

	AtomPtr getMajor()
	{
		return _major.lock();
	}

	AtomPtr getMinor()
	{
		return _minor.lock();
	}

	AtomPtr getHeavyAlign()
	{
		return _heavyAlign.lock();
	}

	AtomPtr getLightAlign()
	{
		return _lightAlign.lock();
	}

	static double getBondLength(void *object)
	{
		return static_cast<Bond *>(object)->_bondLength;
	}

	static void setBondLength(void *object, double length)
	{
		static_cast<Bond *>(object)->_bondLength = length;
	}

	void addDownstreamAtom(AtomPtr influenced);
	void setAlignmentAtoms(AtomPtr heavyAlign, AtomPtr lightAlign);
	virtual FFTPtr getDistribution();
public:
	static double getVoxelValue(void *obj, double x, double y, double z);

private:
	BondGeometryType _minorGeometry;

	AtomWkr _major;
	AtomWkr _minor;

	AtomWkr _heavyAlign;
	AtomWkr _lightAlign;

	double _bondLength;
	double _torsionRadians;

	mat3x3 _torsionBasis;

	std::vector<AtomWkr> _downstreamAtoms;
};

#endif /* defined(__vagabond__Bond__) */
