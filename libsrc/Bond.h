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
#include <string>
#include "Model.h"
#include <vector>
#include "mat3x3.h"
#include "Distributor.h"
#include <iostream>

typedef enum
{
	BondGeometryNone,
	BondGeometryTetrahedral,
} BondGeometryType;

class Bond : public Model
{
public:
	Bond(AtomPtr major, AtomPtr minor);
	void activate(AtomGroupPtr group = AtomGroupPtr());

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

	double getTorsion()
	{
		return _torsionRadians;
	}

	mat3x3 getTorsionBasis()
	{
		return _torsionBasis;
	}

	void setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign);
	virtual FFTPtr getDistribution();
	virtual vec3 getPosition();

	virtual std::string getClassName()
	{
		return "Bond";
	}

	static double getTorsion(void *object)
	{
		return static_cast<Bond *>(object)->_torsionRadians;
	}

	static void setTorsion(void *object, double value)
	{
		static_cast<Bond *>(object)->_torsionRadians = value;
	}

	int downstreamAtomCount()
	{
		return _downstreamAtoms.size();
	}

	AtomPtr downstreamAtom(int i)
	{
		return _downstreamAtoms[i].lock();
	}

	int downstreamAtomNum(AtomPtr atom)
	{
		for (int i = 0; i < _downstreamAtoms.size(); i++)
		{
			if (_downstreamAtoms[i].lock() == atom)
			{
				return i;
			}
		}

		return -1;
	}

	bool isUsingTorsion()
	{
		return _usingTorsion;
	}

	bool isNotJustForHydrogens();
	
protected:
	static double getVoxelValue(void *obj, double x, double y, double z);

private:
	BondGeometryType _minorGeometry;

	AtomWkr _major;
	AtomWkr _minor;

	AtomWkr _heavyAlign;
	AtomWkr _lightAlign;

	double _bondLength;
	double _torsionRadians;

	bool _activated;

	/* Bond direction only used when a torsion angle can't be
	 * calculated because it's connected to an Absolute PDB.
	 * Otherwise use as a reference for torsion matrix updates. */
	vec3 _bondDirection;

	/* Will aim to define torsion basis as:
	   x: along line of 0º torsion angle.
	   y: completes the right-handed coordinate system
	   z: along bond direction, from heavy-to-light alignment atoms.
	 */
	mat3x3 _torsionBasis;
	mat3x3 makeTorsionBasis(vec3 _specificDirection, double *angle = NULL);

	bool _usingTorsion;
};

#endif /* defined(__vagabond__Bond__) */
