//
//  Bond.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include <ostream>
#include <sstream>
#include <iomanip>

#include "Bond.h"
#include "Atom.h"
#include "fftw3d.h"
#include "Shouter.h"
#include "AtomGroup.h"
#include "Element.h"
#include "Absolute.h"
#include "maths.h"
#include "Monomer.h"
#include "Molecule.h"
#include "Anisotropicator.h"
#include "RefinementNelderMead.h"
#include "Options.h"

void Bond::initialize()
{
	_anchored = false;
	_usingTorsion = false;
	_activated = false;
	_dampening = Options::getDampen();
	_bondLength = 0;
	_changedPos = true;
	_changedSamples = true;
	_refineBondAngle = false;
	_refineFlexibility = true;
	_fixed = false;
	_splitBlock = false;
	_occupancy = 1.0;
	_occMult = 1.0;
	_resetOccupancy = false;
	_anisotropyExtent = 0.0;
	_bondDirection = empty_vec3();
	_magicAxis = empty_vec3();
	double initialKick = Options::getKick();
	_kick = initialKick;
	_torsion = 0;
	_phi = 0;
	_phi = 0;
	_geomRatio = 0;
	_circlePortion = -10;
	_expectedAngle = 0;
}

Bond::Bond()
{
	initialize();
}

Bond::Bond(AtomPtr major, AtomPtr minor, int group)
{
	initialize();
	_major = major;
	_minor = minor;

	_disabled = (!major || !minor);

	if (_disabled)
	{
		return;
	}
	
	if (getMinor()->getModel())
	{
		if (getMinor()->getModel()->getClassName() == "Bond")
		{
			std::cout << "Warning!" << std::endl;
		}
	}

	ModelPtr upModel = getParentModel();

	if (upModel->isBond())
	{
		ToBondPtr(upModel)->addDownstreamBond(this, group);
	}

	vec3 majorPos = getMajor()->getInitialPosition();
	vec3 minorPos = getMinor()->getInitialPosition();

	MoleculePtr molecule = getMinor()->getMolecule();
	setMolecule(molecule);

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	_bondLength = vec3_length(difference);

	deriveBondLength();
	vec3_set_length(&difference, _bondLength);

	_bondDirection = difference;
	deriveBondLength();
	deriveBondAngle();
	deriveCirclePortion();
	deriveTorsionAngle();

	if (upModel->getClassName() == "Absolute")
	{
		ToAbsolutePtr(upModel)->addNextAtom(minor);
	}
}

Bond::Bond(Bond &other)
{
	_occupancy = other._occupancy;
	_occMult = other._occMult;
	_usingTorsion = other._usingTorsion;
	_activated = other._activated;
	_major = other._major;
	_minor = other._minor;
	_resetOccupancy = other._resetOccupancy;

	for (size_t i = 0; i < other.downstreamBondGroupCount(); i++)
	{
		BondGroupPtr group = BondGroupPtr(new BondGroup(i));
		_bondGroups.push_back(group);
	}

	_fixed = other._fixed;
	_disabled = other._disabled;
	_anchored = other._anchored;
	_absolute = other._absolute;
	_molecule = other._molecule;

	_refineBondAngle = other._refineBondAngle;
	_refineFlexibility = other._refineFlexibility;
	_dampening = other._dampening;
	_bondLength = other._bondLength;
	_changedPos = true;
	_changedSamples = true;
	_bondDirection = other._bondDirection;
	_heavyAlign = other._heavyAlign;
	_lightAlign = other._lightAlign;
	
	_torsion = other._torsion;
	_kick = other._kick;
	_phi = other._phi;
	_psi = other._psi;
	_circlePortion = other._circlePortion;
	_geomRatio = other._geomRatio;
	_expectedAngle = other._expectedAngle;
}

void Bond::deriveBondLength()
{
	AtomType type1 = getMajor()->getGeomType();
	AtomType type2 = getMinor()->getGeomType();

	GeomTable table = GeomTable::getGeomTable();
	double length = table.getBondLength(type1, type2);

	if (length > 0)
	{
		_bondLength = length;
	}
	else if (getMinor()->getElement()->electronCount() <= 1)
	{
		_bondLength = 0.967;
	}
	else
	{
		std::cout << "Unassigned bond length for " << getMajor()->shortDesc()
		<< " to " << getMinor()->shortDesc() << "." << std::endl;
	}
}

void Bond::setTorsionAngleFrom(AtomPtr one, AtomPtr two, AtomPtr three,
                               AtomPtr four)
{
	vec3 hPos = one->getInitialPosition();
	vec3 maPos = two->getInitialPosition();
	vec3 miPos = three->getInitialPosition();
	vec3 lPos = four->getInitialPosition();

	makeTorsionBasis(hPos, maPos, miPos, lPos,
	                 &_torsion);
}

void Bond::setHeavyAlign(AtomPtr atom, bool from_sister)
{
	_heavyAlign = atom;

	/* if the parent is not available, this is not meant to be
	 * able to generate a torsion angle. Set as heavy alignment atom
	 * for getManyPositions(). */
	if (!getParentModel()->isBond())
	{
		return;
	}
	
	BondPtr parent = ToBondPtr(getParentModel());
	
	if (!parent->getParentModel()->isBond())
	{	
		/* Grandparent is absolute - we can generate a torsion angle
		 * from this heavy atom */

		AtomPtr one = atom;
		AtomPtr two = parent->getMajor();
		AtomPtr three = getMajor();
		AtomPtr four = getMinor();

		setTorsionAngleFrom(one, two, three, four);
		
		/* Rederive circle portion */
		deriveCirclePortion();
	}
	
	if (from_sister)
	{
		return;
	}
	
	/* We may also have sister bonds who should also get this alignment
	 * atom. */

	for (int i = 0; i < parent->downstreamBondGroupCount(); i++)
	{
		for (int j = 0; j < parent->downstreamBondCount(i); j++)
		{
			BondPtr sis = parent->downstreamBond(i, j);

			if (&*sis == this)
			{
				continue;
			}
			
			sis->setHeavyAlign(atom, true);
		}
	}
}

void Bond::deriveTorsionAngle()
{
	if (!getParentModel()->isBond())
	{
		/* It's an absolute - ignore */
		return;
	}
	
	BondPtr parent = ToBondPtr(getParentModel());
	
	if (!parent->getParentModel()->isBond())
	{	
		/* Grandparent is absolute - needs dealing with specially */
		return;
	}
	
	BondPtr grandparent = ToBondPtr(parent->getParentModel());
	
	AtomPtr one = grandparent->getMajor();
	AtomPtr two = parent->getMajor();
	AtomPtr three = getMajor();
	AtomPtr four = getMinor();
	
	setTorsionAngleFrom(one, two, three, four);
}

void Bond::deriveBondAngle()
{
	if (!getParentModel() || !getParentModel()->isBond())
	{
		return;
	}

	BondPtr parent = ToBondPtr(getParentModel());
	AtomPtr pMajor = parent->getMajor();
	
	_expectedAngle = -1;
	double angle = Atom::getAngle(getMinor(), getMajor(), pMajor);
	
	/* In some cases, may not be able to assign, in which case
	 * 	we must fish from the original model */
	if (angle < 0)
	{
		vec3 aPos = pMajor->getInitialPosition();
		vec3 bPos = getMajor()->getInitialPosition();
		vec3 cPos = getMinor()->getInitialPosition();
		vec3 cbDiff = vec3_subtract_vec3(cPos, bPos);
		vec3 abDiff = vec3_subtract_vec3(aPos, bPos);

		/* Derive the angle from the model */
		double angle = vec3_angle_with_vec3(abDiff, cbDiff);
		_expectedAngle = angle;
	}
	
	/* Either way, we need to store the geometry ratio (means
	 *  we don't need to constantly take the tan of something) */
	double ratio = tan(angle - M_PI / 2);
	
	_geomRatio = ratio;
}

double Bond::empiricalCirclePortion(Bond *lastBond)
{
	vec3 hPos;

	BondPtr parent = ToBondPtr(getParentModel());
	if (!parent->getParentModel()->isBond())
	{	
		/* Grandparent is absolute - needs dealing with specially */
		hPos = getHeavyAlign()->getInitialPosition();
	}
	else
	{
		BondPtr grandparent = ToBondPtr(parent->getParentModel());
		hPos = grandparent->getMajor()->getInitialPosition();
	}

	vec3 maPos = parent->getMajor()->getInitialPosition();
	vec3 miPos = getMajor()->getInitialPosition();
	vec3 lPos = getMinor()->getInitialPosition();

	double newAngle = 0;
	makeTorsionBasis(hPos, maPos, miPos, lPos, &newAngle);

	double oldAngle = lastBond->_torsion;
	double increment = newAngle - oldAngle;

	increment /= (2 * M_PI);

	return increment;
}

void Bond::deriveCirclePortion()
{
	if (!getParentModel() || !getParentModel()->isBond())
	{
		return;
	}
	
	BondPtr parent = ToBondPtr(getParentModel());
	
	/* We'll be in the same group as the last group */
	int groups = parent->downstreamBondGroupCount();

	/* First we check to see if there is a sister bond, given that
	 * 	we have not yet been added to the parent */
	int count = parent->downstreamBondCount(groups - 1);
	
	if (getMinor()->getElement()->electronCount() == 1)
	{
		return;
	}
	
	if (count == 1)
	{
		/* We are the first bond. */
		_usingTorsion = true;
		_circlePortion = 0;
		return;
	}
	
	/* First bond should not have a circle portion */
	int num = parent->downstreamBondNum(this, NULL);
	
	if (num == 0)
	{
		return;
	}
	
	/* We can't calculate and compare a circle portion if our grandparent
	 * is an Absolute model and the heavy alignment atom is not set.
	 * This function will be called again when it is set, for now we just
	 * exit. */
	if (!parent->getParentModel()->isBond())
	{	
		if (_heavyAlign.expired())
		{
			return;
		}
	}
	
	/* We have a sister bond and must get a circle portion */
	Bond *lastBond = parent->nakedDownstreamBond(groups - 1, count - 2);

	AtomType central = parent->getMinor()->getGeomType();
	AtomType preceding = parent->getMajor()->getGeomType();
	AtomType lastAtom = lastBond->getMinor()->getGeomType();
	AtomType newDownAtom = getMinor()->getGeomType();

	/* Organise angles to rotate y/z around x */
	/* Which means that angle_a should match x axis */
	GeomTable table = GeomTable::getGeomTable();
	double angle_c = table.getBondAngle(preceding, central, lastAtom);
	double angle_b = table.getBondAngle(preceding, central, newDownAtom);
	double angle_a = table.getBondAngle(lastAtom, central, newDownAtom);

	bool ok = true;

	if (angle_a < 0 || angle_b < 0 || angle_c < 0)
	{
		/* We can't do this - we must derive from the model */
		ok = false;
	}
	
	if (ok)
	{
		mat3x3 bondcell = mat3x3_from_unit_cell(1, 1, 1, rad2deg(angle_a),
		                                        rad2deg(angle_b),
		                                        rad2deg(angle_c));

		vec3 xAxis = mat3x3_axis(bondcell, 0);
		vec3 newAtomAxis = mat3x3_axis(bondcell, 1);
		vec3 lastAtomAxis = mat3x3_axis(bondcell, 2);
		double increment = 0;
		mat3x3_closest_rot_mat(lastAtomAxis, newAtomAxis, xAxis, &increment);
		double theoretical = increment / (2 * M_PI);

		/* This angle will always come out positive, so we must compare
		 * to the original placement as well to maintain chirality */
		double empirical = empiricalCirclePortion(lastBond);
		
		/* Now we clamp both the empirical and the theoretical values to
		 * 	between -0.5 and +0.5 for comparison */
		
		if (empirical > 0.5) empirical -= 1;
		if (empirical < -0.5) empirical += 1;
		
		if (theoretical > 0.5) theoretical -= 1;
		if (theoretical < -0.5) theoretical += 1;
		
		/* there may be some remaining ambiguity around 0.5 here? */

		if ((empirical < 0 && theoretical > 0) || empirical > 0 && theoretical < 0)
		{
			theoretical *= -1;
		}
		
		/* Make the portion positive again if necessary */
		
		if (theoretical < 0) theoretical += 1;

		if (theoretical == theoretical)
		{
			_circlePortion = lastBond->_circlePortion + theoretical;
		}
		else
		{
			/* We can't do this - we must derive from the model */
			ok = false;
		}
	}
	
	if (!ok)
	{
		/* We resort to getting the value from the model.
		 * We calculate the torsion angle for this atom, then subtract
		 * that from the first bond and normalise. */

		double increment = empiricalCirclePortion(lastBond);

		_circlePortion = lastBond->_circlePortion + increment;
	}
}

void Bond::resetBondAngles()
{
	deriveBondAngle();
}

void Bond::addDownstreamBond(Bond *bond, int group)
{
	if (_bondGroups.size() <= group)
	{
		for (int i = 0; i < group - _bondGroups.size() + 1; i++)
		{
			int index = _bondGroups.size();
			BondGroupPtr group = BondGroupPtr(new BondGroup(index));
			_bondGroups.push_back(group);
		}
	}

	_bondGroups[group]->addBond(bond);
}

void Bond::setMinor(AtomPtr newMinor)
{
	_minor = newMinor;

	if (_bondLength < 0)
	{
		deriveBondLength();
	}

	vec3 majorPos = getMajor()->getInitialPosition();
	vec3 minorPos = getMinor()->getInitialPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	vec3_set_length(&difference, _bondLength);

	_bondDirection = difference;
}

void Bond::activate()
{
	if (_disabled)
	return;

	getMinor()->setModel(shared_from_this());

	_activated = true;

}

mat3x3 Bond::makeTorsionBasis(vec3 hPos, vec3 maPos,
                              vec3 miPos, vec3 lPos, double *newAngle)
{
	vec3 a2p = vec3_subtract_vec3(hPos, maPos);
	vec3 a2l = vec3_subtract_vec3(lPos, maPos);
	vec3 a2b = vec3_subtract_vec3(miPos, maPos);

	double sql = vec3_sqlength(a2b);
	double dotp = vec3_dot_vec3(a2p, a2b);
	double dotl = vec3_dot_vec3(a2l, a2b);
	double distp = dotp / sql;
	double distl = dotl / sql;

	vec3 a2bp = a2b;
	vec3 a2bl = a2b;
	vec3_mult(&a2bp, distp);
	vec3_mult(&a2bl, distl);
	vec3 heavy_join = vec3_add_vec3(maPos, a2bp);
	vec3 light_join = vec3_add_vec3(maPos, a2bl);

	vec3 reverse_bond = vec3_subtract_vec3(miPos, maPos);
	vec3 xNew = vec3_subtract_vec3(hPos, heavy_join);
	vec3 light = vec3_subtract_vec3(lPos, light_join);
	mat3x3 basis = mat3x3_rhbasis(xNew, reverse_bond);

	if (newAngle)
	{
		double angle = vec3_angle_with_vec3(light, xNew);

		vec3 recross = vec3_cross_vec3(light, xNew);
		double cosine = vec3_cosine_with_vec3(recross, reverse_bond);

		if (cosine > 0)
		{
			angle *= -1;
		}

		if (angle != angle)
		{
			shout_at_helen("Torsion angle is nan for " + shortDesc() + "!");
		}


		*newAngle = angle;
	}

	return basis;
}

FFTPtr Bond::makeDistribution()
{
	return makeRealSpaceDistribution();
}

vec3 Bond::positionFromTorsion(mat3x3 torsionBasis, double angle,
                               double ratio, vec3 start)
{
	double myLength = getBondLength(this);     // my bond length!

	/* Calculate the right ratio of x-to-z from the major atom. */
	vec3 atomWrtMajor = make_vec3(1, 0, ratio);

	/* Set these to adhere to the correct bond length. */
	vec3_set_length(&atomWrtMajor, myLength);

	/* Rotate the bond so it lines up with the torsion angle. */
	vec3 zAxis = make_vec3(0, 0, 1);
	mat3x3 torsion_turn = mat3x3_unit_vec_rotation(zAxis, -angle);
	mat3x3_mult_vec(torsion_turn, &atomWrtMajor);

	/* Reset the basis vectors so that they line up with the previous
	* bond and the torsion angle is correct. */
	mat3x3_mult_vec(torsionBasis, &atomWrtMajor);

	/* Add this to the major atom position! */
	vec3 final = vec3_add_vec3(start, atomWrtMajor);

	return final;
}

mat3x3 Bond::getMagicMat(vec3 direction)
{
	vec3_set_length(&direction, 1.);
	mat3x3 rot = make_mat3x3();

	if (_phi != 0 || _psi != 0)
	{
		rot = mat3x3_rot_from_angles(_phi, _psi);
	}

	vec3 xAxis = make_vec3(1, 0, 0);
	vec3 zAxis = make_vec3(0, 0, 1);
	
	vec3 cross = vec3_cross_vec3(zAxis, direction);
	vec3_set_length(&cross, 1.);

	/* Find the twizzle to put z axis onto the magic axis (around the x) */
	mat3x3 firstTwizzle = mat3x3_closest_rot_mat(direction, zAxis, cross,
	                                             NULL, true);
	mat3x3 multed = mat3x3_mult_mat3x3(rot, firstTwizzle);

	return multed;
}

typedef struct
{
	mat3x3 basis;
	vec3 curr;
	vec3 next;	
} BondCache;

void Bond::correctTorsionAngles(std::vector<BondSample> *prevs)
{
	const vec3 none = make_vec3(0, 0, 0);

	mat3x3 aveBasis = make_mat3x3();
	vec3 aveStart = make_vec3(0, 0, 0);

	/* This loop gets average positions for the basis of the bond
	 * and the average start position */
	for (size_t i = 0; i < prevs->size(); i++)
	{
		mat3x3_add_mat3x3(&aveBasis, prevs->at(i).basis);
		aveStart = vec3_add_vec3(aveStart, prevs->at(i).start);
	}

	double samples = prevs->size();
	mat3x3_mult_scalar(&aveBasis, 1 / samples);
	vec3_mult(&aveStart, 1 / samples);

	vec3 aveNext = mat3x3_axis(aveBasis, 0);
	
	vec3 crossDir = mat3x3_axis(aveBasis, 1);
	mat3x3 magicMat = getMagicMat(crossDir);
	_magicAxis = mat3x3_axis(magicMat, 2); 
	
	/* Track overall change in order to readjust torsion
	 * at the end */
	double averageModulation = 0;
	
	/* Assume torsion of 0, as real torsion added later */
	for (size_t i = 0; i < prevs->size(); i++)
	{
		mat3x3 thisBasis = (*prevs)[i].basis;
		vec3 prevHeavyPos = (*prevs)[i].old_start;
		vec3 thisPos = prevs->at(i).start;

		/* Difference between perfect and deviant position of major atom */
		vec3 thisDeviation = vec3_subtract_vec3(thisPos, aveStart);
		vec3_set_length(&thisDeviation, 1.);
		
		/* Find out what this deviation is if beam axis is set to z */
		mat3x3_mult_vec(magicMat, &thisDeviation);

		double notZ = sqrt(thisDeviation.y * thisDeviation.y +
		                   thisDeviation.x * thisDeviation.x);
		double tanX = thisDeviation.z / notZ;
		double angle = atan(tanX);
		
		if (thisDeviation.z < 0) angle *= -1;
		
		double dampValue = sin(angle);
		if (dampValue != dampValue)
		{
			dampValue = 0;
		}
		
		double kickValue = cos(angle);
		if (kickValue != kickValue)
		{
			kickValue = 0;
		}

		/* (dampening) We want to correct if the deviation is close to 
		 * the magic angle */
		double rotAngle = 0;
		
		/** Average direction of THIS bond (torsion-ignorant) */
		vec3 thisDir = mat3x3_axis(thisBasis, 2);
		
		/** Actual direction of THIS bond (torsion-ignorant) */
		vec3 thisNext = mat3x3_axis(thisBasis, 0);

		/* Find the best angle for dampening */
		mat3x3_closest_rot_mat(thisNext, aveNext, thisDir, &rotAngle, true);

		if (rotAngle != rotAngle)
		{
			rotAngle = 0;
		}

		double undoBlur = 0;
		undoBlur = rotAngle;
		undoBlur *= fabs(dampValue);
		undoBlur *= fabs(_dampening);

		/* This will only apply for a kicked bond */
		double addBlur = _kick;
		addBlur *= kickValue;

		if (isFixed())
		{
			undoBlur = 0; addBlur = 0;
		}

		double totalBlur = undoBlur + addBlur;
		prevs->at(i).torsion = totalBlur;	
		
		averageModulation += totalBlur;
	}

	averageModulation /= (double)samples;

	for (size_t i = 0; i < prevs->size(); i++)
	{
		prevs->at(i).torsion -= averageModulation;
	}

}

double Bond::getBaseTorsion()
{
	ModelPtr model = getParentModel();
	BondPtr prevBond = boost::static_pointer_cast<Bond>(model);
	int myGroup = -1;
	double torsionNumber = prevBond->downstreamBondNum(this,
	                                                   &myGroup);

	BondPtr parent = ToBondPtr(getParentModel());
	BondPtr sisBond = parent->downstreamBond(myGroup, 0);
	/* May be myself */
	double baseTorsion = sisBond->_torsion;
	return baseTorsion;
}

std::vector<BondSample> *Bond::getManyPositions()
{
	std::vector<BondSample> *newSamples;

	newSamples = &_storedSamples;

	if (!_changedSamples)
	{
		return &_storedSamples;
	}

	ModelPtr model = getParentModel();

	newSamples->clear();

	if (model->getClassName() == "Absolute")
	{
		std::vector<BondSample> *absPos = model->getManyPositions();
		mat3x3 magicMat = getMagicMat(_bondDirection);
		
		_magicAxis = mat3x3_axis(magicMat, 2); 

		/* We must be connected to something else, oh well */
		/* Torsion basis must be the same. */

		for (size_t i = 0; i < absPos->size(); i++)
		{
			vec3 majorPos = (*absPos)[i].start;
			vec3 heavyPos = getHeavyAlign()->getInitialPosition();
			vec3 none = {0, 0, 1};
			vec3 actualMajor = getMajor()->getAbsolutePosition();
			vec3 start = vec3_subtract_vec3(majorPos, _bondDirection);
			vec3 perfectStart = vec3_subtract_vec3(actualMajor, _bondDirection);

			mat3x3 newBasis = makeTorsionBasis(heavyPos, actualMajor, perfectStart, none);

			vec3 majorDev = vec3_subtract_vec3(majorPos, actualMajor);
			mat3x3_mult_vec(magicMat, &majorDev);
			double notY = sqrt(majorDev.x * majorDev.x +
			                   majorDev.z * majorDev.z);
			double tanY = majorDev.y / notY;
			double diffValue = sin(atan(tanY));

			if (diffValue != diffValue)
			{
				diffValue = 0;
			}

			BondSample newSample;
			newSample.basis = newBasis;
			newSample.torsion = 0;
			newSample.start = start;
			newSample.old_start = majorPos;
			newSample.occupancy = (*absPos)[i].occupancy;
			newSamples->push_back(newSample);

		}

		return newSamples;
	}

	BondPtr prevBond = boost::static_pointer_cast<Bond>(model);
	int myGroup = -1;
	double torsionNumber = prevBond->downstreamBondNum(this,
	                                                   &myGroup);

	bool nextBondExists = false;
	if (downstreamBondGroupCount())
	{
		nextBondExists = true;
	}

	double totalAtoms = prevBond->downstreamBondCount(myGroup);

	std::vector<BondSample> *prevSamples = prevBond->getManyPositions();
	
	double baseTorsion = getBaseTorsion();

	/* This is just to get a set of angles, no bases */
	double circlePortion = 0;
	double circleAdd = 0;

	circlePortion = _circlePortion;

	if (circlePortion < -9) /* Likely a hydrogen */
	{
		circleAdd += deg2rad(360) * torsionNumber / totalAtoms;
	}
	else
	{
		circleAdd += deg2rad(360) * circlePortion;
	}

	double ratio = _geomRatio;

	std::vector<BondSample> myTorsions;

	bool usingKick = (isUsingTorsion()
	                          && nextBondExists && !isFixed());

	if (usingKick)
	{
		correctTorsionAngles(prevSamples);
	}

	double occTotal = 0;
	newSamples->reserve(prevSamples->size());

	for (size_t i = 0; i < prevSamples->size(); i++)
	{
		double currentTorsion = baseTorsion + circleAdd;
		/* Deviation from correction */
		currentTorsion += prevSamples->at(i).torsion;

		vec3 prevHeavyPos = (*prevSamples)[i].old_start;
		vec3 prevMinorPos = (*prevSamples)[i].start;
		mat3x3 oldBasis = (*prevSamples)[i].basis;

		vec3 myCurrentPos = positionFromTorsion(oldBasis, currentTorsion,
		                                        ratio, prevMinorPos);

		/* Prepping for bending */
		const vec3 none = {0, 0, 1};

		/* New basis for the next bond */
		mat3x3 newBasis = makeTorsionBasis(prevHeavyPos, prevMinorPos,
		                                   myCurrentPos, none);

		BondSample nextSample;
		nextSample.basis = newBasis;
		nextSample.start = myCurrentPos;
		nextSample.old_start = prevMinorPos;
		nextSample.torsion = 0;
		nextSample.occupancy = (_occupancy * prevSamples->at(i).occupancy);
		
		occTotal += nextSample.occupancy;

		newSamples->push_back(nextSample);
	}

	if (_resetOccupancy)
	{
		for (size_t i = 0; i < newSamples->size(); i++)
		{
			newSamples->at(i).occupancy /= occTotal;
		}
	}

	_changedSamples = false;

	return newSamples;
}

bool Bond::isNotJustForHydrogens()
{
	if (getMinor()->getElement()->electronCount() > 1)
	{
		return true;
	}

	if (downstreamBondGroupCount() == 0)
	{
		return false;
	}

	for (size_t i = 0; i < downstreamBondCount(0); i++)
	{
		AtomPtr atom = downstreamAtom(i, 0);
		if (atom->getElement()->electronCount() > 1)
		{
			return true;
		}
	}

	return false;
}

void Bond::propagateChange(int depth, bool refresh)
{
	if (_anchored)
	{
		return;
	}

	/* Will force recalculation of final positions */
	Model::propagateChange(depth, refresh);

	/* Iterative now */
	std::vector<BondPtr> propagateBonds;
	propagateBonds.push_back(ToBondPtr(shared_from_this()));
	int count = 0;

	for (size_t k = 0; k < propagateBonds.size(); k++)
	{
		BondPtr bond = propagateBonds[k];

		for (size_t j = 0; j < bond->downstreamBondGroupCount(); j++)
		{
			for (size_t i = 0; i < bond->downstreamBondCount(j); i++)
			{
				AtomPtr atom = bond->downstreamAtom(j, i);
				ModelPtr model = atom->getModel();
				BondPtr bond = ToBondPtr(model);

				propagateBonds.push_back(bond);
			}
		}

		bond->_changedPos = true;
		bond->_changedSamples = true;
		bond->_recalcDist = true;
		bond->Model::propagateChange(depth, refresh);

		if (depth >= 0 && count > depth)
		{
			break;
		}

		count++;
	}

	if (!refresh)
	{
		return;
	}

	for (size_t k = 0; k < propagateBonds.size(); k++)
	{
		BondPtr bond = propagateBonds[k];
		{
			bond->getFinalPositions();
		}
	}
}

ModelPtr Bond::getParentModel()
{
	AtomPtr atom = getMajor();
	ModelPtr model = atom->getModel();

	return model;
}

void Bond::setCirclePortion(void *object, double value)
{
	Bond *bond = static_cast<Bond *>(object);
	bond->_circlePortion = value / (2 * M_PI);
}

double Bond::getCirclePortion(void *object)
{
	Bond *bond = static_cast<Bond *>(object);
	return bond->_circlePortion * (2 * M_PI);
}

bool Bond::connectsAtom(std::string name)
{
	return (getMinor()->getAtomName() == name ||
	        getMajor()->getAtomName() == name);

}

double Bond::getBendAngle(void *object)
{
	Bond *bond = static_cast<Bond *>(object);
	double angle = atan(bond->_geomRatio) + M_PI / 2;
	return angle;
}

double Bond::getExpectedAngle()
{
	return _expectedAngle;
}

void Bond::setBendAngle(void *object, double value)
{
	Bond *bond = static_cast<Bond *>(object);
	double ratio = tan(value - M_PI / 2);
	bond->_geomRatio = ratio;

	bond->propagateChange(10);
}

BondGroupPtr Bond::bondGroupForBond()
{
	ModelPtr model = getParentModel();
	
	if (!model->isBond())
	{
		return NULL;
	}

	BondPtr parent = ToBondPtr(model);
	int group = -1;
	int num = parent->downstreamBondNum(this, &group);
	
	if (group < 0)
	{
		shout_at_helen("Group less than zero - floating bond?");
	}

	return parent->_bondGroups[group];
}

bool Bond::splitBond()
{
	BondPtr me = ToBondPtr(shared_from_this());
	BondPtr parent = ToBondPtr(getParentModel());
	int last = parent->downstreamBondGroupCount();
	
	int num = parent->downstreamBondNum(this, NULL);
	
	BondPtr dupl = me->duplicateDownstream(parent, last);

	double torsion = getBaseTorsion();
	
	if (num > 0)
	{
		torsion += getCirclePortion(&*me);
		setCirclePortion(&*dupl, 0.);
	}
	
	setTorsion(&*dupl, torsion);
	dupl->setUsingTorsion(true);

	_occupancy /= 2;
	dupl->_occupancy /= 2;


	propagateChange();

	return true;
}

void Bond::copyParamsFromFirstGroup(BondPtr copyFrom, int groupNum)
{
	/* Set the torsion angle to be the same as the parent */
	_torsion = copyFrom->_torsion;
	_circlePortion = copyFrom->_circlePortion;
	_kick = copyFrom->_kick;
}

BondPtr Bond::duplicateDownstream(BondPtr newParent, int groupNum)
{
	/* duplBond will be a child of newParent */
	BondPtr duplBond = BondPtr(new Bond(*this));

	/* Dealing with atom business */
	AtomPtr duplAtom = AtomPtr();
	
	/* We look for other atoms of the same name which haven't been tied up */
	std::string search = getMinor()->getAtomName();
	AtomList list = getMinor()->getMonomer()->findAtoms(search);

	/* Do we have something in the list which has an Absolute model?
	* if not, we will create a new one. */
	for (size_t i = 0; i < list.size(); i++)
	{
		if (list[i].expired()) continue;
		AtomPtr atom = list[i].lock();
		ModelPtr model = atom->getModel();

		if (model->isAbsolute())
		{
			duplAtom = atom;
			break;
		}
	}

	/* In the case where we have not found an existing atom...*/
	if (!duplAtom)
	{
		duplAtom = AtomPtr(new Atom(*getMinor()));
		/* Lower case */
		char conformer[] = "a";

		for (size_t i = 0; i < list.size(); i++)
		{
			conformer[0]++;
		}

		duplAtom->setAlternativeConformer(conformer);
		duplAtom->setFromPDB(false);
		getMinor()->getMolecule()->addAtom(duplAtom);
	}
	/* Dealt with appropriate atom business */

	/* Connect this duplicated minor atom to the polymer and monomer.
	 * The molecule should already have happened in both cases above */
	duplAtom->inheritParents();
	duplAtom->setModel(duplBond);
	
	/* Set the duplicate bond's major atom to the duplicated parent's 
	 * minor atom */
	duplBond->setMajor(newParent->getMinor());
	duplBond->setMinor(duplAtom);
	
	/* Add the new minor atom to the correct group in the parent bond */
	newParent->addDownstreamBond(&*duplBond, groupNum);
	
	/* If this is the end of the chain, stop now */
	if (!downstreamBondGroupCount())
	{
		return duplBond;
	}
	
	/* Need to spread the duplication to the entire next set of bonds
	 * (currently no support for duplicating already-branched bonds!! FIXME */
	for (size_t i = 0; i < downstreamBondCount(0); i++)
	{
		if (_splitBlock)
		{
			_splitBlock = false;
			/* Need to return occupancy to full here on out */
			ModelPtr model = downstreamAtom(0, i)->getModel();
			ToBondPtr(model)->_resetOccupancy = true;
			
			continue;
		}

		AtomPtr nextAtom = downstreamAtom(0, i);
		ModelPtr nextModel = nextAtom->getModel();

		if (!nextModel->isBond())
		{
			continue;
		}

		BondPtr nextBond = ToBondPtr(nextAtom->getModel());
		nextBond->duplicateDownstream(duplBond, 0);
	}

	return duplBond;
}

void Bond::setOccupancy(void *object, double value)
{
	Bond *bond = static_cast<Bond *>(object);
	bond->_occupancy = value;

	static_cast<Bond *>(object)->propagateChange();
}

std::string Bond::shortDesc()
{
	std::ostringstream stream;
	stream << getMajor()->getMonomer()->getResidueNum() <<
	getMajor()->getAtomName() << "-" << getMinor()->getAtomName();
	if (getMinor()->getAlternativeConformer().length())
	{
		stream << "_" + getMinor()->getAlternativeConformer();
	}

	_shortDesc = stream.str();
	return stream.str();
}

std::string Bond::description()
{
	std::ostringstream stream;
	stream << "Bond: " << shortDesc() << std::endl;
	stream << "Bond length: " << _bondLength << " Å" << std::endl;
	stream << "Bond torsion angle: "
	<< _torsion << std::endl;
	stream << "Bond downstream groups: ("
	<< downstreamBondGroupCount() << "):" << std::endl;
	stream << "Bond downstream atoms (first) ("
	<< downstreamBondCount(0) << "):" << std::endl;

	for (size_t i = 0; i < downstreamBondGroupCount(); i++)
	{
		for (size_t j = 0; j < downstreamBondCount(i); j++)
		{
			stream << "\t" << downstreamAtom(i, j)->shortDesc()
			<< "(" << &*downstreamAtom(i, j) << ")" << std::endl;
		}
	}

	stream << "Major: " << getMajor()->shortDesc() << "("
	<< &*getMajor() << ")" << std::endl;
	stream << "Minor: " << getMinor()->shortDesc() << "("
	<< &*getMinor() << ")" << std::endl;



	return stream.str();
}

double Bond::getMeanSquareDeviation()
{
	getAnisotropy(true);
	return _isotropicAverage * 8 * M_PI * M_PI;
}

void Bond::resetBondDirection()
{
	vec3 majorPos = getMajor()->getAbsolutePosition();
	vec3 minorPos = getMinor()->getAbsolutePosition();

	if (getParentModel()->isBond())
	{
		ToBondPtr(getParentModel())->getDistribution();
		majorPos = ToBondPtr(getParentModel())->getAbsolutePosition();
		getDistribution();
		minorPos = getAbsolutePosition();
	}

	_bondDirection = vec3_subtract_vec3(majorPos, minorPos);
	vec3_set_length(&_bondDirection, _bondLength);
}

bool Bond::isRefinable()
{
	return isNotJustForHydrogens() && !isFixed();
}

bool Bond::test()
{
	bool ok = true;

	/* Test of geometry for multiple downstream atoms */
	for (size_t i = 0; i < downstreamBondGroupCount(); i++)
	{
		for (size_t j = -1; j < downstreamBondCount(i); j++)
		{
			AtomPtr atom1 = getMajor();
			if (j >= 0)
			{
				atom1 = downstreamAtom(i, j);
			}

			if (ToBondPtr(atom1->getModel())->getRefineBondAngle())
			{
				continue;
			}

			for (size_t k = 0; k < downstreamBondCount(i); k++)
			{
				AtomPtr atom3 = downstreamAtom(i, (k + 1) % downstreamBondCount(i));

				if (ToBondPtr(atom3->getModel())->getRefineBondAngle())
				{
					continue;
				}

				double angle = Atom::getAngle(atom1, getMinor(), atom3);

				if (angle < 0)
				{
					/* Has no geometric entry */
					continue;
				}

				/* Comparing conformers which should not be matched! */
				if (atom1->getAlternativeConformer() != atom3->getAlternativeConformer())
				{
					continue;
				}
				//
				vec3 pos1 = atom1->getModel()->getManyPositions()->at(0).start;
				vec3 pos2 = getMinor()->getModel()->getManyPositions()->at(0).start;
				vec3 pos3 = atom3->getModel()->getManyPositions()->at(0).start;

				double realAngle = vec3_angle_from_three_points(pos1, pos2, pos3);

				if (angle > deg2rad(90))
				{
					angle -= deg2rad(90);
				}

				if (realAngle > deg2rad(90))
				{
					realAngle -= deg2rad(90);
				}

				double diff = fabs(realAngle - angle);

				if (diff > 1e-3)
				{
					std::cout << shortDesc() << " angle "
					<< atom1->getAtomName() << "-" << getMinor()->getAtomName() << "-"
					<< atom3->getAtomName() << "\t"
					<< 90 + rad2deg(realAngle) << "\t" << 90 + rad2deg(angle) << std::endl;
				}

				ok *= (diff < 1e-6);
			}
		}
	}

	return ok;
}

void Bond::recalculateTorsion(AtomPtr heavy, double value)
{
	if (!_heavyAlign.expired() && _heavyAlign.lock() != heavy)
	{
		heavy->getModel()->getFinalPositions();
		_heavyAlign.lock()->getModel()->getFinalPositions();
		vec3 newHPos = heavy->getModel()->getAbsolutePosition();
		vec3 origHPos = _heavyAlign.lock()->getModel()->getAbsolutePosition();
		getMinor()->getModel()->getFinalPositions();
		getMajor()->getModel()->getFinalPositions();
		vec3 miPos = getMinor()->getAbsolutePosition();
		vec3 maPos = getMajor()->getAbsolutePosition();

		double unwantedTorsion = 0;
		makeTorsionBasis(origHPos, maPos, miPos, newHPos, &unwantedTorsion);

		value += unwantedTorsion;
	}

	setTorsion(this, value);
}


double Bond::getEffectiveOccupancy()
{
	std::vector<BondSample> *samples = getManyPositions();
	double total = 0;

	for (size_t i = 0; i < samples->size(); i++)
	{
		total += samples->at(i).occupancy;
	}

	return total;
}

void Bond::addProperties()
{
	addReference("minor", getMinor());
	addReference("major", getMajor());
	addReference("heavy", getHeavyAlign());
	addReference("light", getLightAlign());

	addDoubleProperty("length", &_bondLength);    
	addDoubleProperty("torsion", &_torsion);    
	addDoubleProperty("exp_angle", &_expectedAngle);    
	addDoubleProperty("circle", &_circlePortion);    
	addDoubleProperty("ratio", &_geomRatio);    
	addDoubleProperty("kick", &_kick);    
	addDoubleProperty("phi", &_phi);    
	addDoubleProperty("psi", &_psi);    
	addDoubleProperty("dampening", &_dampening);
	addDoubleProperty("occupancy", &_occupancy);
	addDoubleProperty("occ_mult", &_occMult);
	addBoolProperty("occ_reset", &_resetOccupancy);

	addBoolProperty("fixed", &_fixed);
	addBoolProperty("anchored", &_anchored);
	addBoolProperty("using_torsion", &_usingTorsion);
	addBoolProperty("refine_bond_angle", &_refineBondAngle);
	addBoolProperty("refine_flexibility", &_refineFlexibility);
	addBoolProperty("activated", &_activated);
	addBoolProperty("disabled", &_disabled);

	addVec3Property("bond_direction", &_bondDirection);
	addVec3Property("magic_axis", &_magicAxis);
	
	for (int i = 0; i < downstreamBondGroupCount(); i++)
	{
		addChild("bond_group", _bondGroups[i]);
	}
	
	for (int i = 0; i < extraTorsionSampleCount(); i++)
	{
		addReference("extra_sample", extraTorsionSample(i));
	}
	
	Model::addProperties();
}

void Bond::addObject(ParserPtr object, std::string category)
{
	if (category == "bond_group")
	{
		BondGroupPtr group = ToBondGroupPtr(object);
		_bondGroups.push_back(group);	
	} 
}

void Bond::linkReference(ParserPtr object, std::string category)
{
	AtomPtr atom = ToAtomPtr(object);

	if (category == "minor")
	{
		_minor = atom;
	}
	else if (category == "major")
	{
		_major = atom;
	}
	else if (category == "heavy")
	{
		_heavyAlign = atom;
	}
	else if (category == "light")
	{
		_lightAlign = atom;
	}
	else if (category == "extra_sample")
	{
		addExtraTorsionSample(atom);
	} 
}

void Bond::postParseTidy()
{	

}
