//
//  Bond.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

/**
* \class Bond
*
* \brief The Bond class is the most complex type of Model, which controls a
* single Atom. Bonds create a network and recursively inherit and build upon
* one another.
*
* These consist of several linked atoms. A bond controls a single atom, the minor atom (available through setMinor() and getMinor()). This immediately branches from the major atom (setMajor() and getMajor()) which is controlled by the parent bond. Groups of downstream atoms (one for every conformer) can be accessed using downstreamAtom(int, int).
*
* A small number of bonds will have a parent Absolute model (those next to the anchor atom of a Polymer) and therefore getParentModel() should not be assumed to be a Bond.
*/


#ifndef __vagabond__Bond__
#define __vagabond__Bond__

#include <stdio.h>
#include <string>
#include <vector>
#include "mat3x3.h"
#include "Distributor.h"
#include <iostream>
#include "Atom.h"
#include "BondGroup.h"
#include "Sampler.h"
#include "Model.h"
#include <mutex>
#include "charmanip.h"

#define INITIAL_KICK 0.01
#define INITIAL_DAMPENING 0.08

class Anisotropicator;

class Bond : public Model
{
public:
	Bond(AtomPtr major, AtomPtr minor, int group = 0);
	Bond(Bond &other);
	Bond();
	virtual ~Bond() {};

	/**
	* After a bond has been initialised with a major and minor (pre-existing)
	* atom, it must be activated. It must be activated after all parent bonds
	* have been activated. It then sets up the references to the other atoms.
	* */
	void activate();
	ModelPtr reverse(BondPtr upstreamBond);
	void reverseDownstreamAtoms(int group);
	void resetBondDirection();
	void setTorsionAngleFrom(AtomPtr one, AtomPtr two, AtomPtr three,
	                         AtomPtr four);
	
	BondPtr shared_from_this()
	{
		return ToBondPtr(Model::shared_from_this());
	}
	
	/** Returns true if the torsion angle should be refined. */
	bool isRefinable();
	
	/** Returns true if the bond angles and lengths for the first atom agree
	* with default geometry or are otherwise refinable to non-default values.
	* */
	bool test();
	double getEffectiveOccupancy();
	
	/** Resets bond angles to the default from geometry (n.b. used for fixing
	* disulphides. */
	void resetBondAngles();

	/**
	*  Major atom from which bond is drawn. Bond does not dictate position
	*  of the major atom.
	*/
	AtomPtr getMajor()
	{
		return _major.lock();
	}

	/**
	*  Bond directly controls the distribution of the minor atom.
	*/
	AtomPtr getMinor()
	{
		return _minor.lock();
	}

	void setMajor(AtomPtr newMajor)
	{
		_major = newMajor;
	}

	/** Sets minor atom and re-derives bond length for the bond. */
	void setMinor(AtomPtr newMinor);

	/**
	*  This is the first (preceding) atom on which the torsion angle is
	*  calculated.
	*  \return heavy alignment atom.	
	*/
	AtomPtr getHeavyAlign()
	{
		return _heavyAlign.lock();
	}
	
	void setHeavyAlign(AtomPtr atom, bool from_sister = false);

	/**
	*  This is the second atom on which the torsion angle is
	*  calculated. Perfect alignment of heavy atom and light atom will
	*  be a torsion angle of 0.
	*  \return light alignment atom.	
	*/
	AtomPtr getLightAlign()
	{
		return _lightAlign.lock();
	}
	
	/**
	* 	If an atom of a given name is part of this bond.
	* \return If the minor or major atom contains this atom name (e.g. "CA").*/
	bool connectsAtom(std::string testName);

	static double getBondLength(void *object)
	{
		return static_cast<Bond *>(object)->_bondLength;
	}

	static void setBondLength(void *object, double length)
	{
		static_cast<Bond *>(object)->_bondLength = length;
	}

	virtual FFTPtr makeDistribution();
	
	virtual AtomPtr getAtom()
	{
		return getMinor();
	}

	virtual std::string getClassName()
	{
		return "Bond";
	}

	static void setKick(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_kick = value;
		static_cast<Bond *>(object)->propagateChange(16);
	}

	static double getKick(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_kick;
	}

	static double getTorsion(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_torsion;
	}
	
	static void setTorsion(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_torsion = value;
		static_cast<Bond *>(object)->propagateChange(16);
	}

	static double getMagicPsi(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_psi;
	}

	static void setMagicPsi(void *object, double angle)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_psi = angle;
		static_cast<Bond *>(object)->propagateChange(16);
	}

	static double getMagicPhi(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_phi;
	}

	static void setMagicPhi(void *object, double angle)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_phi = angle;
		static_cast<Bond *>(object)->propagateChange(16);
	}

	static double getOccupancy(void *object)
	{
		Bond *bond = static_cast<Bond *>(object);
		return bond->_occupancy;
	}

	static void setOccupancy(void *object, double value);

	static double getDampening(void *object)
	{
		return static_cast<Bond *>(object)->_dampening;
	}

	static void setDampening(void *object, double value)
	{
		Bond *bond = static_cast<Bond *>(object);
		bond->_dampening = value;
		bond->propagateChange(16);
	}

	/**
	* 	Returns the bond angle for a bond, between upstream->major->minor.
	*   \param object C pointer to the bond object, cast as void *.
	* 	\return angle in radians.
	*/
	static double getBendAngle(void *object);
	
	/**
	* 	Sets the bond angle for a bond between upstream->major->minor.
	* 	\param object C pointer to the bond object, cast as void *.
	* 	\param value value of new angle in radians.
	*/
	static void setBendAngle(void *object, double value);

	/**
	*	Adds a downstream atom to a given group. Groups represent split
	* 	downstream bonds.
	*   \param atom to add to downstream
	*   \param group to add atom to in the current bond
	*   \param skipGeometry if set to true, geometry is not recalculated
	*   	from default values.
	*/
	void addDownstreamAtom(AtomPtr atom, int group, bool skipGeometry = false);

	/**
	* 	\return Number of downstream conformers.
	*/
	size_t downstreamBondGroupCount()
	{
		return _bondGroups.size();
	}

	/**
	*  Returns a minor atom from the immediate descendent bond.	
	*  \param group group to fetch atom from.
	*  \param i number of the atom in the array within the group.
	*  \return pointer to downstream atom.
	*/
	AtomPtr downstreamAtom(int group, int i)
	{
		return downstreamBond(group, i)->getMinor();
	}

	/**
	*   \param group which group to query.
	* 	\return Number of downstream atoms in a given group (conformer).
	*/
	size_t downstreamBondCount(int group)
	{
		return _bondGroups[group]->bondCount();
	}


	/**
	* 	Returns the appropriate group/number for a given bond, expected
	* 	to be called by a direct descendant of a given bond.
	*	\param bond One of the bonds in one of the groups of this bond.
	*	\param group pointer to be filled with the expected group.
	*	\return -1 if bond isn't found, and the bond number if found.
	*		
	*/
	int downstreamBondNum(Bond *down, int *group)
	{
		for (int j = 0; j < downstreamBondGroupCount(); j++)
		{
			for (int i = 0; i < downstreamBondCount(j); i++)
			{
				if (nakedDownstreamBond(j, i) == down)		
				{
					if (group)
					{
						*group = j;
					}

					return i;
				}
			}
		}

		if (group)
		{
			*group = -1;
		}

		return -1;
	}

	void setUsingTorsion(bool use)
	{
		_usingTorsion = use;
	}

	bool isUsingTorsion()
	{
		return _usingTorsion;
	}


	virtual bool hasExplicitPositions()
	{
		return true;
	}

	/** 
	* 	Bond torsion does something other than change hydrogen atom position.
	*/
	bool isNotJustForHydrogens();

	/**
	*  This is updated when the geometry is recalculated. This is the expected
	*  angle from default geometry.
	*  \return expected bond angle in radians between upstream->major->minor.
	*/
	double getExpectedAngle();
	
	/**
	*  Geometry ratio for direct calculation of atom position.
	*/
	double getGeomRatio()
	{
		return _geomRatio;
	}

	
	/**
	*  In cases where there are multiple downstream bonds, this returns the
	*  offset of the torsion from the first atom in the downstream atom array.
	*  \return number of radians to add to the first torsion angle.
	*/
	double getCirclePortion(int n, int i)
	{
		return _circlePortion;
	}

	static double getCirclePortion(void *object);
	static void setCirclePortion(void *object, double value);

	void setGeomRatio(double value)
	{
		_geomRatio = value;
	}

	BondGroupPtr getBondGroup(int i)
	{
		return _bondGroups[i];
	}

	void setOccupancyMult(double mult)
	{
		_occMult = mult;
		propagateChange();
	}

	std::string description();
	std::string shortDesc();
	std::string getPDBContribution();
	ModelPtr getParentModel();

	/**
	* 	 Splits all downstream atoms and creates new copies of atoms and
	* 	 bonds. Sets initial torsion angle difference to 180 degrees of first
	* 	 bond, and divides the occupancies by 2 by default.
	* 	 \param start if set, downstream atoms will only be duplicated from
	* 	 this position in the array.
	*	\return Always returns true, at the moment.
	*/
	bool splitBond();

	void setFixed(bool fixed)
	{
		_fixed = fixed;
	}

	bool isFixed()
	{
		return _fixed;
	}

	void addExtraTorsionSample(AtomPtr atom)
	{
		_extraTorsionSamples.push_back(atom);
	}

	int extraTorsionSampleCount()
	{
		return _extraTorsionSamples.size();
	}

	/** Returns the B factor (function is a misnomer). */
	virtual double getMeanSquareDeviation();

	AtomPtr extraTorsionSample(int i)
	{
		return _extraTorsionSamples[i].lock();
	}


	/** Will define torsion basis as:
	* x: along line of 0º torsion angle.
	* y: completes the right-handed coordinate system
	* z: along bond direction, from heavy-to-light alignment atoms.
	*/
	mat3x3 makeTorsionBasis(vec3 hPos, vec3 maPos,
	                        vec3 miPos, vec3 lPos, double *newAngle = NULL);

	/** Can determine a new torsion angle with a different heavy atom.
	* 	Would be useful in cases where the chain needs to be reversed */
	void recalculateTorsion(AtomPtr heavy, double value);

	virtual void propagateChange(int depth = -1, bool refresh = false);
	std::vector<BondSample> *getManyPositions();

	static void useMutex()
	{
		_useMutex = true;
	}

	double getMultOccupancy()
	{
		if (_resetOccupancy) return 1;
		return _occupancy * _occMult;
	}
	
	void setRefineFlexibility(bool value = true)
	{
		_refineFlexibility = value;
	}
	
	bool getRefineFlexibility()
	{
		return _refineFlexibility && !isFixed();
	}

	void setRefineBondAngle(bool value = true)
	{
		_refineBondAngle = value;
	}

	bool getRefineBondAngle()
	{
		return _refineBondAngle && !isFixed();
	}

	static void encodeBondGroup(void *bond, void *bondGroup,
	                            std::ostream &stream, int indent);
	static char *decodeBondGroup(void *bond, void *bondGroup, char *block);

	void setSplitBlock(int block = 1)
	{
		_splitBlock = block;
	}
protected:

	AtomWkr _minor;

	virtual std::string getParserIdentifier()
	{
		return "bond_" + shortDesc();
	}

	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category);
	virtual void addObject(ParserPtr object, std::string category);
	virtual void postParseTidy();    
	friend class StateValue;

private:
	void initialize();
	double getBaseTorsion();

	std::string _shortDesc;

	AtomWkr _major;

	AtomWkr _heavyAlign;
	AtomWkr _lightAlign;

	double _bondLength;

	/* Downstream groups of bonds */
	std::vector<BondGroupPtr> _bondGroups;
	
	
	Bond *nakedDownstreamBond(int group, int i)
	{
		return _bondGroups[group]->bond(i);
	}
	
	BondPtr downstreamBond(int group, int i)
	{
		return nakedDownstreamBond(group, i)->shared_from_this();
	}

	/* Returns upstream bond group pertaining to this bond. */
	BondGroupPtr bondGroupForBond();

	std::vector<AtomWkr> _extraTorsionSamples;
	std::vector<BondSample> _storedSamples;

	/* Dampening should be associated with a bond group - woops */
	double _dampening;
	bool _activated;
	double _occupancy;
	double _occMult;
	double _torsion;
	double _kick;
	double _phi;
	double _psi;
	double _circlePortion;
	double _geomRatio;
	double _expectedAngle;
	
	bool _resetOccupancy;
	
	/* If blocked, do not duplicate downstream */
	int _splitBlock;

	/* Should not be refined */
	bool _fixed;

	/* Had a non-NULL atom input as major or minor */
	bool _disabled;

	/* Has been set as an anchor, will not respond to 'propagate change'*/
	bool _anchored;

	/* Grab bond length from the atom types of major/minor */
	void deriveBondLength();
	void deriveBondAngle();
	void deriveCirclePortion();
	void deriveTorsionAngle();
	double empiricalCirclePortion(Bond *lastBond);

	void addDownstreamBond(Bond *bond, int group);

	/* Bond direction only used when a torsion angle can't be
	* calculated because it's connected to an Absolute PDB.
	* Otherwise use as a reference for torsion matrix updates. */
	vec3 _bondDirection;
	
	vec3 _magicAxis;

	/** Supply deviations of correct torsion angles into prevs->torsion */
	void correctTorsionAngles(std::vector<BondSample> *prevs);

	vec3 positionFromTorsion(mat3x3 torsionBasis, double angle,
	                         double ratio, vec3 start);

	void copyParamsFromFirstGroup(BondPtr copyFrom, int groupNum);
	BondPtr duplicateDownstream(BondPtr newBranch, int groupNum);
	bool _usingTorsion;

	/* Flag to say whether recalculation should occur */
	bool _changedPos, _changedSamples;
	bool _refineBondAngle;
	bool _refineFlexibility;

	mat3x3 getMagicMat(vec3 direction);


};

#endif /* defined(__vagabond__Bond__) */
