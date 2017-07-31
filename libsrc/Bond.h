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

typedef struct
{
	mat3x3 basis;
	vec3 start;
	double torsion;
	double occupancy;
} BondSample;

class Bond : public Model
{
public:
	Bond(AtomPtr major, AtomPtr minor);
	void activate(AtomGroupPtr group = AtomGroupPtr(),
				  AbsolutePtr inherit = AbsolutePtr());

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
	virtual vec3 getStaticPosition();
	std::vector<BondSample> getManyPositions(bool staticAtom = false);
	virtual std::string getClassName()
	{
		return "Bond";
	}

	static void setTorsionBlur(void *object, double value)
	{
		static_cast<Bond *>(object)->_torsionBlur = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getTorsionBlur(void *object)
	{
		return static_cast<Bond *>(object)->_torsionBlur;
	}

	static double getTorsion(void *object)
	{
		return static_cast<Bond *>(object)->_torsionRadians;
	}

	static void setTorsion(void *object, double value)
	{
		static_cast<Bond *>(object)->_torsionRadians = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getTorsionNextBlur(void *object)
	{
		return static_cast<Bond *>(object)->_torsionBlurFromPrev;
	}

	static void setTorsionNextBlur(void *object, double value)
	{
		static_cast<Bond *>(object)->_torsionBlurFromPrev = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getBendBlur(void *object)
	{
		return static_cast<Bond *>(object)->_bendBlur;
	}

	static void setBendBlur(void *object, double value)
	{
		static_cast<Bond *>(object)->_bendBlur = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getBendAngle(void *object);
	static void setBendAngle(void *object, double value);

	void setBendTowards(AtomWkr bendToAtom)
	{
		_bendToAtom = bendToAtom;
	}

	void addDownstreamAtom(AtomPtr atom);

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

	double getGeomRatio(int i)
	{
		return _downRatios[i];
	}

	void setGeomRatio(int i, double value)
	{
		_downRatios[i] = value;
	}

	void setAbsoluteInheritance(AbsolutePtr abs)
	{
		_absInherit = abs;
	}

	std::string description();
	std::string getPDBContribution();
	ModelPtr getParentModel();

protected:
	static double getVoxelValue(void *obj, double x, double y, double z);

private:
	AtomWkr _major;
	AtomWkr _minor;

	AtomWkr _heavyAlign;
	AtomWkr _lightAlign;
	AtomWkr _bendToAtom;

	double _bondLength;
	double _torsionRadians;
	double _torsionBlur;
	double _torsionBlurFromPrev;
	double _bendBlur;

	bool _activated;

	/* Bond direction only used when a torsion angle can't be
	 * calculated because it's connected to an Absolute PDB.
	 * Otherwise use as a reference for torsion matrix updates. */
	vec3 _bondDirection;

	std::vector<double> _downRatios;

	/* Will aim to define torsion basis as:
	   x: along line of 0º torsion angle.
	   y: completes the right-handed coordinate system
	   z: along bond direction, from heavy-to-light alignment atoms.
	 */
	mat3x3 _torsionBasis;
	mat3x3 makeTorsionBasis(vec3 hPos, vec3 maPos,
							vec3 miPos, vec3 lPos, double *newAngle = NULL);

	vec3 positionFromTorsion(mat3x3 torsionBasis, double angle,
							 double ratio, vec3 start);
	std::vector<BondSample> sampleMyAngles(double angle, double sigma,
										   double interval = ANGLE_SAMPLING);
	std::vector<BondSample> getCorrelatedAngles(BondSample prev,
												double lastTorsion,
												double angle, double blur);

	void propagateChange();
	bool _usingTorsion;

	/* Flag to say whether recalculation should occur */
	bool _changedPos, _changedSamples;
	vec3 _lastPosition;
	std::vector<BondSample> _lastSamples;

	AbsolutePtr _absInherit;
};

#endif /* defined(__vagabond__Bond__) */
