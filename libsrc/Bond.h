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
#include "Atom.h"

typedef struct
{
	mat3x3 basis;
	vec3 start;
	vec3 old_start;
	double torsion;
	double occupancy;
} BondSample;

typedef enum
{
	BondSampleThorough,
	BondSampleStatic,
	BondSampleMonteCarlo
} BondSampleStyle;

typedef struct
{
	AtomWkr atom;
	double geomRatio;
	double circlePortion;
} AtomValue;

typedef struct
{
	std::vector<AtomValue> atoms;
	double torsionAngle;
	double torsionBlur;
	double occupancy;
	std::vector<BondSample> storedSamples;
	std::vector<BondSample> staticSample;
	std::vector<BondSample> singleStateSample;
	std::vector<AtomWkr> extraTorsionSamples;
	bool _changedSamples;
} BondGroup;

class Bond : public Model
{
public:
	Bond(AtomPtr major, AtomPtr minor, int group = 0);
	Bond(Bond &other);
	void activate(AtomGroupPtr group = AtomGroupPtr(),
				  AtomPtr inherit = AtomPtr());

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
		return _bondGroups[group].torsionAngle;
	}

	void setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign);
	virtual FFTPtr getDistribution();
	virtual vec3 getStaticPosition();
	std::vector<BondSample> *getManyPositions(BondSampleStyle style);
	virtual std::string getClassName()
	{
		return "Bond";
	}

	static void setTorsionBlur(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_bondGroups[bond->_activeGroup].torsionBlur = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getTorsionBlur(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_bondGroups[bond->_activeGroup].torsionBlur;
	}

	static double getTorsion(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_bondGroups[bond->_activeGroup].torsionAngle;
	}

	static void setTorsion(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_bondGroups[bond->_activeGroup].torsionAngle = value;
		static_cast<Bond *>(object)->propagateChange();
	}

	static double getOccupancy(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_bondGroups[bond->_activeGroup].occupancy;
	}

	static void setOccupancy(void *object, double value);
	
	static double getTorsionNextBlur(void *object)
	{
		return static_cast<Bond *>(object)->_torsionBlurFromPrev;
	}

	static void setTorsionNextBlur(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_torsionBlurFromPrev = std::max(-1., -fabs(value));
		bond->propagateChange();
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
		return _bondGroups[group].atoms.size();
	}

	int downstreamAtomGroupCount()
	{
		return _bondGroups.size();
	}

	AtomPtr downstreamAtom(int group, int i)
	{
		return _bondGroups[group].atoms[i].atom.lock();
	}

	int downstreamAtomNum(AtomPtr atom, int *group)
	{
		for (int j = 0; j < downstreamAtomGroupCount(); j++)
		{
			for (int i = 0; i < downstreamAtomCount(j); i++)
			{
				if (downstreamAtom(j, i) == atom)
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
		return _bondGroups[n].atoms[i].geomRatio;
	}

	double getCirclePortion(int n, int i)
	{
		return _bondGroups[n].atoms[i].circlePortion;
	}

	void setGeomRatio(int n, int i, double value)
	{
		_bondGroups[n].atoms[i].geomRatio = value;
	}

	void setAbsoluteInheritance(AtomPtr abs)
	{
		_absInherit = abs;
	}

	void setActiveGroup(int newGroup)
	{
		_activeGroup = newGroup;
		propagateChange();
	}

	int getActiveGroup()
	{
		return _activeGroup;
	}

	std::string description();
	std::string shortDesc();
	std::string getPDBContribution();
	ModelPtr getParentModel();
	bool splitBond();

	void setFixed(bool fixed)
	{
		_fixed = fixed;
	}

	bool isFixed()
	{
		return _fixed;
	}

	void addExtraTorsionSample(AtomPtr atom, int group)
	{
		_bondGroups[group].extraTorsionSamples.push_back(atom);
	}

	int extraTorsionSampleCount(int group)
	{
		return _bondGroups[group].extraTorsionSamples.size();
	}

	AtomPtr extraTorsionSample(int group, int i)
	{
		return _bondGroups[group].extraTorsionSamples[i].lock();
	}

protected:

private:
	AtomWkr _major;
	AtomWkr _minor;

	AtomWkr _heavyAlign;
	AtomWkr _lightAlign;
	AtomWkr _bendToAtom;

	double _bondLength;

	/* Bond groups */
	std::vector<BondGroup> _bondGroups;
	std::vector<AtomWkr> _extraTorsionSamples;

	double _torsionBlurFromPrev;
	double _bendBlur;
	bool _activated;
	int _activeGroup;
	bool _fixed;

	/* Grab bond length from the atom types of major/minor */
	void deriveBondLength();
	/* And a given bond angle for downstream atom n */
	void deriveBondAngle(int group, int n);

	/* Bond direction only used when a torsion angle can't be
	 * calculated because it's connected to an Absolute PDB.
	 * Otherwise use as a reference for torsion matrix updates. */
	vec3 _bondDirection;
	mat3x3 _torsionBasis;

	/* Will aim to define torsion basis as:
	   x: along line of 0º torsion angle.
	   y: completes the right-handed coordinate system
	   z: along bond direction, from heavy-to-light alignment atoms.
	 */
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
	virtual void propagateChange();
	bool _usingTorsion;

	/* Flag to say whether recalculation should occur */
	bool _changedPos, _changedSamples;

	AbsolutePtr getAbsInheritance();

	AtomPtr _absInherit;
	FFTPtr _fftAbs;
};

#endif /* defined(__vagabond__Bond__) */
