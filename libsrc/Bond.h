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

typedef struct
{
	AtomWkr atom;
	double geomRatio;
} AtomValue;

typedef struct
{
	std::vector<AtomValue> atoms;
	double torsionAngle;
	double torsionBlur;
	double occupancy;
} BondGroup;

typedef std::vector<AtomWkr> AtomList;

class Bond : public Model
{
public:
	Bond(AtomPtr major, AtomPtr minor, int group = 0);
	Bond(Bond &other);
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

	void setMajor(AtomPtr newMajor)
	{
		_major = newMajor;
	}

	void setMinor(AtomPtr newMinor);

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

	double getTorsion(int group)
	{
		return _torsionAngles[group];
	}

	mat3x3 getTorsionBasis()
	{
		return _torsionBasis;
	}


	void setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign);
	virtual FFTPtr getDistribution();
	virtual vec3 getStaticPosition();
	std::vector<BondSample> getManyPositions(bool staticAtom = false,
											 bool singleState = false,
											 int torsionGroup = 0);
	virtual std::string getClassName()
	{
		return "Bond";
	}

	static void setTorsionBlur(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_torsionBlurs[bond->_activeGroup] = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getTorsionBlur(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_torsionBlurs[bond->_activeGroup];
	}

	static double getTorsion(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_torsionAngles[bond->_activeGroup];
	}

	static void setTorsion(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_torsionAngles[bond->_activeGroup] = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getTorsionNextBlur(void *object)
	{
		return static_cast<Bond *>(object)->_torsionBlurFromPrev;
	}

	static void setTorsionNextBlur(void *object, double value)
	{
		static_cast<Bond *>(object)->_torsionBlurFromPrev = std::max(-1., -fabs(value));
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

	void addDownstreamAtom(AtomPtr atom, int group);

	int downstreamAtomCount(int group)
	{
		return _downstreamAtoms[group].size();
	}

	int downstreamAtomGroupCount()
	{
		return _downstreamAtoms.size();
	}

	AtomPtr downstreamAtom(int group, int i)
	{
		return _downstreamAtoms[group][i].lock();
	}

	int downstreamAtomNum(AtomPtr atom, int *group)
	{
		for (int j = 0; j < _downstreamAtoms.size(); j++)
		{
			for (int i = 0; i < _downstreamAtoms[j].size(); i++)
			{
				if (_downstreamAtoms[j][i].lock() == atom)
				{
					*group = j;
					return i;
				}
			}
		}

		*group = -1;
		return -1;
	}

	bool isUsingTorsion()
	{
		return _usingTorsion;
	}

	bool isNotJustForHydrogens();

	double getGeomRatio(int n, int i)
	{
		return _downRatios[n][i];
	}

	void setGeomRatio(int n, int i, double value)
	{
		_downRatios[n][i] = value;
	}

	void setAbsoluteInheritance(AbsolutePtr abs)
	{
		_absInherit = abs;
	}

	void setActiveGroup(int newGroup)
	{
		_activeGroup = newGroup;
	}

	std::string description();
	std::string getPDBContribution();
	ModelPtr getParentModel();
	bool splitBond();

protected:
	static double getVoxelValue(void *obj, double x, double y, double z);

private:
	AtomWkr _major;
	AtomWkr _minor;

	std::vector<AtomList> _downstreamAtoms;

	AtomWkr _heavyAlign;
	AtomWkr _lightAlign;
	AtomWkr _bendToAtom;

	double _bondLength;

	/* Bond groups */
	std::vector<BondGroup> bondGroups;

	/* To be "bond-grouped" */
	std::vector<double> _torsionAngles, _torsionBlurs;
	std::vector<std::vector<double> > _downRatios;
	double _torsionBlurFromPrev;
	double _bendBlur;
	/* END: To be "bond-grouped" */

	bool _activated;
	int _activeGroup;

	/* Grab bond length from the atom types of major/minor */
	void deriveBondLength();
	/* And a given bond angle for downstream atom n */
	void deriveBondAngle(int group, int n);

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
	mat3x3 makeTorsionBasis(vec3 hPos, vec3 maPos,
							vec3 miPos, vec3 lPos, double *newAngle = NULL);

	vec3 positionFromTorsion(mat3x3 torsionBasis, double angle,
							 double ratio, vec3 start);
	std::vector<BondSample> sampleMyAngles(double angle, double sigma,
										   bool singleState = false);
	std::vector<BondSample> getCorrelatedAngles(BondSample prev,
												double lastTorsion,
												double angle, double blur,
												bool singleState = false,
												int group = 0);

	void duplicateDownstream(BondPtr newBranch, int groupNum);
	void propagateChange();
	bool _usingTorsion;

	/* Flag to say whether recalculation should occur */
	bool _changedPos, _changedSamples;
	vec3 _lastPosition;
	std::vector<BondSample> _lastSamples;

	AbsolutePtr _absInherit;
};

#endif /* defined(__vagabond__Bond__) */
