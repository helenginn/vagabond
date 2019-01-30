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
#include "Whack.h"
#include "Twist.h"
#include "Anchor.h"
#include "maths.h"
#include "Monomer.h"
#include "Polymer.h"
#include "Options.h"

Bond::~Bond()
{

}

void Bond::initialize()
{
	_usingTorsion = false;
	_bondLength = 0;
	_changedSamples = true;
	_split = false;
	_refineBondAngle = false;
	_refineFlexibility = true;
	_fixed = false;
	_splitBlock = false;
	_occupancy = 1.0;
	_resetOccupancy = false;
	_anisotropyExtent = 0.0;
	_magicAxis = empty_vec3();
	double initialKick = Options::getKick();
	_kick = initialKick;
	_torsion = 0;
	_phi = 3.14/2;
	_psi = 0;
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

	deriveBondLength();
	deriveBondAngle();
	deriveCirclePortion();
	deriveTorsionAngle();
	
}

Bond::Bond(Bond &other)
{
	_occupancy = other._occupancy;
	_usingTorsion = other._usingTorsion;
	_major = other._major;
	_minor = other._minor;
	_resetOccupancy = false;
	_magicAxis = other._magicAxis;

	for (size_t i = 0; i < other.downstreamBondGroupCount(); i++)
	{
		BondGroupPtr group = BondGroupPtr(new BondGroup(i));
		_bondGroups.push_back(group);
	}

	_fixed = other._fixed;
	_disabled = other._disabled;
	_absolute = other._absolute;
	_molecule = other._molecule;
	_splitBlock = false;

	_refineBondAngle = other._refineBondAngle;
	_refineFlexibility = other._refineFlexibility;
	_bondLength = other._bondLength;
	_changedSamples = true;
	_heavyAlign = other._heavyAlign;
	_split = other._split;
	
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
	else if (getMinor()->getElectronCount() <= 1)
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
	vec3 hPos = empty_vec3();
	if (one)
	{
		hPos = one->getInitialPosition();
	}

	vec3 maPos = two->getInitialPosition();
	vec3 miPos = three->getInitialPosition();
	vec3 lPos = four->getInitialPosition();

	mat3x3 basis = makeTorsionBasis(hPos, maPos, miPos, lPos,
	                 &_torsion);
}

void Bond::setHeavyAlign(AtomPtr atom, bool from_sister)
{
	_heavyAlign = atom;

	if (!getParentModel()->isBond())
	{
		deriveTorsionAngle();
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
	if (getParentModel()->isAnchor())
	{
		AnchorPtr anchor = ToAnchorPtr(getParentModel());
		/* It's an absolute - make relative to origin, to have
		 * something which can be reproduced */
		AtomPtr one = getHeavyAlign();
		AtomPtr two = anchor->getOtherAtom(getMinor());
		AtomPtr three = getMajor();
		AtomPtr four = getMinor();

		setTorsionAngleFrom(one, two, three, four);
		return;
	}
	
	BondPtr parent = ToBondPtr(getParentModel());
	
	if (!parent->getParentModel()->isBond())
	{	
		/* Grandparent is absolute - needs dealing with specially,
		 * but can be relative. */
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
	if (!getParentModel())
	{
		return;
	}
	
	AtomPtr pMajor;
	
	if (getParentModel()->isAnchor())
	{
		AnchorPtr parent = ToAnchorPtr(getParentModel());
		pMajor = parent->getOtherAtom(getMinor());
	}
	else
	{
		BondPtr parent = ToBondPtr(getParentModel());
		pMajor = parent->getMajor();
	}
	
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
		if (getParentModel()->isAnchor())
		{
			setUsingTorsion(true);
			_circlePortion = 0;
		}

		return;
	}
	
	BondPtr parent = ToBondPtr(getParentModel());
	
	/* We'll be in the same group as the last group */
	int groups = parent->downstreamBondGroupCount();

	/* First we check to see if there is a sister bond, given that
	 * 	we have not yet been added to the parent */
	int count = parent->downstreamBondCount(groups - 1);
	
	if (getMinor()->getElectronCount() == 1)
	{
		return;
	}
	
	if (count == 1)
	{
		/* We are the first bond. */
		setUsingTorsion(true);
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

}

void Bond::activate()
{
	if (_disabled)
	return;

	getMinor()->setModel(shared_from_this());
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

mat3x3 Bond::getMagicMat(mat3x3 basis)
{
	mat3x3 rot = make_mat3x3();

	if (_phi != 0 || _psi != 0)
	{
		rot = mat3x3_rot_from_angles(_phi, _psi);
	}

	mat3x3 multed = mat3x3_mult_mat3x3(rot, basis);

	return mat3x3_transpose(multed);
}

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
	
	mat3x3 magicMat = getMagicMat(aveBasis);
	_magicAxis = mat3x3_axis(magicMat, 2); 
//	magicMat = make_mat3x3();
	
	/* Track overall change in order to readjust torsion
	 * at the end */
	double averageModulation = 0;
	
	/* Assume torsion of 0, as real torsion added later */
	for (size_t i = 0; i < prevs->size(); i++)
	{
		mat3x3 thisBasis = (*prevs)[i].basis;
		vec3 thisPos = prevs->at(i).start;

		/* Difference between perfect and deviant position of major atom */
		vec3 thisDeviation = vec3_subtract_vec3(thisPos, aveStart);
		vec3_set_length(&thisDeviation, 1.);
		
		/* Find out what this deviation is if sensitive axis is set to z */
		mat3x3_mult_vec(magicMat, &thisDeviation);

		double notZ = (thisDeviation.y * thisDeviation.y +
		               thisDeviation.x * thisDeviation.x);
		double c = thisDeviation.z * thisDeviation.z / notZ;

		double sinAlpha = sqrt(c / (c + 1));
		
//		if (thisDeviation.z < 0) sinAlpha *= -1;

		double kickValue = sinAlpha;
		
		/* Just the notZ */
		if (kickValue != kickValue)
		{
			kickValue = 0;
		}

		/* Baseline kick multiplied by kickValue */
		double addBlur = _kick * kickValue;

		if (isFixed())
		{
			addBlur = 0;
		}

		double totalBlur = addBlur;
		prevs->at(i).torsion = totalBlur;	
		prevs->at(i).kickValue = kickValue;
		
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
	ExplicitModelPtr model = getParentModel();
	int myGroup = -1;
	double torsionNumber = model->downstreamBondNum(this, &myGroup);

	/* May be myself */
	BondPtr sisBond = shared_from_this();

	if (torsionNumber > 0)
	{
		BondPtr parent = ToBondPtr(getParentModel());
		sisBond = parent->downstreamBond(myGroup, 0);
	}

	double baseTorsion = sisBond->_torsion;
	
	if (torsionNumber == 0 && model->isBond())
	{
		BondPtr parent = ToBondPtr(getParentModel());
		
		if (parent->hasTwist())
		{
			double angle = parent->getTwist()->getTorsionCorrection();
			baseTorsion += angle;
		}
	}
	
	return baseTorsion;
}

std::vector<BondSample> *Bond::getManyPositions(void *)
{
	if (!_changedSamples)
	{
		return &_storedSamples;
	}

	ExplicitModelPtr model = getParentModel();

	std::vector<ExplicitModelPtr> parents;

	/* Go back through the parental models until we find one
	 * which does not have _changedSamples == true */
	while (true)
	{
		if (model->changedSamples())
		{
			parents.push_back(model);
		}
		else
		{
			break;
		}
		
		if (model->isAnchor())
		{
			break;
		}
		else if (model->isBond())
		{
			model = ToBondPtr(model)->getParentModel();
		}
	}

	/* Now we have a list of parents in reverse order; call each
	 * individually to get refreshed positions. */

	for (int i = parents.size() - 1; i >= 0; i--)
	{
		/* If the parent is actually an anchor, we skip it; it will
		 * get appropriately called by an actual bond */
		if (parents[i]->isBond())
		{
			ToBondPtr(parents[i])->getManyPositionsPrivate();
		}
	}
	
	/* Now we can call it on ourselves */
	
	std::vector<BondSample> *samples = getManyPositionsPrivate();
	
	/* And return it */
	
	return samples;
}
	
std::vector<BondSample> *Bond::getManyPositionsPrivate()
{
	if (!_changedSamples)
	{
		return &_storedSamples;
	}

	ExplicitModelPtr model = getParentModel();
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

	_storedSamples = *prevBond->getManyPositions(&*getMinor());
	
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
	                          && !isFixed());

	if (usingKick)
	{
		correctTorsionAngles(&_storedSamples);
	}

	double occTotal = 0;
	/* vec3 to be ignored */
	const vec3 none = {0, 0, 1};

	for (size_t i = 0; i < _storedSamples.size(); i++)
	{
		double currentTorsion = baseTorsion + circleAdd;
		/* Deviation from correction */
		currentTorsion += _storedSamples[i].torsion;

		vec3 prevHeavyPos = _storedSamples[i].old_start;
		vec3 prevMinorPos = _storedSamples[i].start;
		mat3x3 oldBasis = _storedSamples[i].basis;

		vec3 myCurrentPos = positionFromTorsion(oldBasis, currentTorsion,
		                                        ratio, prevMinorPos);

		/* New basis for the next bond */
		mat3x3 newBasis = makeTorsionBasis(prevHeavyPos, prevMinorPos,
		                                   myCurrentPos, none);

		_storedSamples[i].basis = newBasis;
		_storedSamples[i].start = myCurrentPos;
		_storedSamples[i].old_start = prevMinorPos;
		_storedSamples[i].occupancy *= _occupancy;
		
		occTotal += _storedSamples[i].occupancy;
	}

	if (_resetOccupancy)
	{
		for (size_t i = 0; i < _storedSamples.size(); i++)
		{
			_storedSamples[i].occupancy /= occTotal;
		}
	}

	_changedSamples = false;
	
	sanityCheck();

	return &_storedSamples;
}

bool Bond::isNotJustForHydrogens()
{
	if (getMinor()->getElectronCount() > 1)
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
		if (atom->getElectronCount() > 1)
		{
			return true;
		}
	}

	return false;
}

AtomGroupPtr Bond::makeAtomGroup(BondPtr endBond)
{
	AtomGroupPtr group = AtomGroupPtr(new AtomGroup());

	std::vector<BondPtr> propagateBonds;
	propagateBonds.push_back(shared_from_this());

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
				
				if (bond != endBond)
				{
					propagateBonds.push_back(bond);
				}

			}
		}
		
		AtomPtr minor = bond->getMinor();
		group->addAtom(minor);
	}
	
	return group;
}

void Bond::propagateChange(int depth, bool refresh)
{
	/* Will force recalculation of final positions */
	Model::propagateChange(depth, refresh);

	/* Iterative now */
	std::vector<BondInt> propagateBonds;
	BondInt pair;
	pair.bond = shared_from_this();
	pair.num = depth;
	propagateBonds.push_back(pair);
	
	/* Propagate down sisters */
	ModelPtr pModel = getParentModel();
	
	if (pModel && pModel->isBond())
	{
		BondPtr parent = ToBondPtr(pModel);

		int group = -1;
		int num = parent->downstreamBondNum(this, &group);
		
		if (group >= 0 && num == 0)
		{
			for (int j = 1; j < parent->downstreamBondCount(group); j++)
			{
				BondPtr downBond = parent->downstreamBond(group, j);
				
				BondInt pair;
				pair.bond = downBond;
				pair.num = depth;
				propagateBonds.push_back(pair);
			}
		}
	}

	for (size_t k = 0; k < propagateBonds.size(); k++)
	{
		BondPtr bond = propagateBonds[k].bond;
		int curr_depth = propagateBonds[k].num;

		if ((curr_depth >= 0 && depth >= 0) || depth < 0)
		{
			for (size_t j = 0; j < bond->downstreamBondGroupCount(); j++)
			{
				for (size_t i = 0; i < bond->downstreamBondCount(j); i++)
				{
					BondPtr newBond = bond->downstreamBond(j, i);
					BondInt pair;
					pair.bond = newBond;
					pair.num = curr_depth - 1;

					propagateBonds.push_back(pair);
				}
			}
		}

		bond->_changedSamples = true;
		bond->_recalcDist = true;
		/* not recursive call - this fixes other problems */
		bond->Model::propagateChange(depth, refresh);
	}

	if (!refresh)
	{
		return;
	}

	for (size_t k = 0; k < propagateBonds.size(); k++)
	{
		BondPtr bond = propagateBonds[k].bond;
		{
			bond->getFinalPositions();
		}
	}
}

ExplicitModelPtr Bond::getParentModel()
{
	AtomPtr atom = getMajor();
	ExplicitModelPtr model = atom->getExplicitModel();

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
		return BondGroupPtr();
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

void Bond::sanityCheck()
{
	if (_torsion != _torsion)
	{
		std::cout << "Torsion angle is nan" << std::endl;
		std::cout << description() << std::endl;
	}
	
	ExplicitModel::sanityCheck();
}

void Bond::destroy(bool start)
{
	for (int i = 0; i < downstreamBondGroupCount(); i++)
	{
		for (int j = 0; j < downstreamBondCount(i); j++)
		{
			downstreamBond(i, j)->destroy(false);
		}
	}
	
	if (start)
	{
		BondPtr parent = ToBondPtr(getParentModel());
		int group = -1;
		int num = parent->downstreamBondNum(this, &group);

		parent->_bondGroups[group]->removeBond(num);
	}
	
	PolymerPtr polymer = getMinor()->getMonomer()->getPolymer();
	polymer->removeAtom(getMinor());
	
	if (getParentModel()->isAnchor())
	{
		warn_user("Pruning too far... reached anchor point.");
	}
}

void Bond::checkForSplits(AtomGroupPtr polymer)
{
	AtomPtr minor = getMinor();
	AtomList list = polymer->findAtoms(minor->getAtomName(), 
	                                   minor->getResidueNum());
	
	if (list.size() > 1)
	{
		BondPtr success = BondPtr();
		for (int i = 0; i < list.size() - 1; i++)
		{
			BondPtr result = splitBond(true);
			
			if (result)
			{
				success = result;
			}
		}

		if (success)
		{
			for (int i = 0; i < list.size(); i++)
			{
				AtomPtr atom = list[i];
				if (atom->getModel()->isBond())
				{
					BondPtr bond = ToBondPtr(atom->getModel());
					double origOcc = atom->getOriginalOccupancy();
					Bond::setOccupancy(&*bond, origOcc);
				}
			}
		}

	}
	
	for (size_t j = 0; j < downstreamBondGroupCount(); j++)
	{
		for (size_t i = 0; i < downstreamBondCount(j); i++)
		{
			BondPtr bond = downstreamBond(j, i);
			bond->checkForSplits(polymer);	
		}
	}
}

BondPtr Bond::splitBond(bool onlyExisting)
{
	BondPtr me = ToBondPtr(shared_from_this());
	BondPtr parent = ToBondPtr(getParentModel());
	int last = parent->downstreamBondGroupCount();
	
	int num = parent->downstreamBondNum(this, NULL);
	
	BondPtr dupl = me->duplicateDownstream(parent, last, onlyExisting);
	
	if (dupl)
	{
		if (getMinor()->getAlternativeConformer().length() == 0)
		{
			getMinor()->setAlternativeConformer("a");
		}
	}
	
	if (!dupl && onlyExisting)
	{
		return BondPtr();
	}

	double torsion = getBaseTorsion();
	
	if (num > 0)
	{
		torsion += getCirclePortion(&*me);
		setCirclePortion(&*dupl, 0.);
		dupl->setFixed(true);
	}
	
	setTorsion(&*dupl, torsion);
	dupl->setUsingTorsion(true);

	_occupancy /= 2;
	dupl->_occupancy /= 2;

	dupl->setSplit(true);
	setSplit(true);

	propagateChange();

	return dupl;
}

void Bond::copyParamsFromFirstGroup(BondPtr copyFrom, int groupNum)
{
	/* Set the torsion angle to be the same as the parent */
	_torsion = copyFrom->_torsion;
	_circlePortion = copyFrom->_circlePortion;
	_kick = copyFrom->_kick;
}

BondPtr Bond::duplicateDownstream(BondPtr newParent, int groupNum,
                                  bool onlyExisting)
{
	/* Dealing with atom business */
	AtomPtr duplAtom = AtomPtr();
	
	/* We look for other atoms of the same name which haven't been tied up */
	std::string search = getMinor()->getAtomName();
	AtomList list = getMinor()->getMonomer()->findAtoms(search);

	/* Do we have something in the list which has an Absolute model?
	* if not, we will create a new one. */
	for (size_t i = 0; i < list.size(); i++)
	{
		AtomPtr atom = list[i];
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
		if (onlyExisting)
		{
			return BondPtr();
		}

		CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
		duplAtom = AtomPtr(new Atom(*getMinor()));
		/* Lower case */
		char conformer[] = "a";

		if (getMinor()->getAlternativeConformer().length() == 0)
		{
			getMinor()->setAlternativeConformer(conformer);
		}
		
		for (size_t i = 0; i < list.size(); i++)
		{
			conformer[0]++;
		}

		duplAtom->setAlternativeConformer(conformer);
		duplAtom->setFromPDB(false);
		int num = crystal->issueAtomNumber();
		duplAtom->setAtomNum(num);
		getMinor()->getMolecule()->addAtom(duplAtom);
	}

	/* Dealt with appropriate atom business */

	/* duplBond will be a child of newParent */
	BondPtr duplBond = BondPtr(new Bond(*this));

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
		if (onlyExisting)
		{
			duplBond->deriveTorsionAngle();
		}

		return duplBond;
	}
	
	/* Need to spread the duplication to the entire next set of bonds
	 * (currently no support for duplicating already-branched bonds!! FIXME */
	for (size_t i = 0; i < downstreamBondCount(0); i++)
	{
		if (_splitBlock)
		{
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
		BondPtr result = nextBond->duplicateDownstream(duplBond, 
		                                               0, onlyExisting);
		
		if (onlyExisting && !result)
		{
			nextBond->_resetOccupancy = true;
		}
	}

	_splitBlock = false;
	
	if (onlyExisting)
	{
		duplBond->deriveTorsionAngle();
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
	
	if (downstreamBondGroupCount())
	{
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
	}

	stream << "Major: " << getMajor()->shortDesc() << "("
	<< &*getMajor() << ")" << std::endl;
	stream << "Minor: " << getMinor()->shortDesc() << "("
	<< &*getMinor() << ")" << std::endl;



	return stream.str();
}

bool Bond::isRefinable()
{
	return isNotJustForHydrogens() && !isFixed();
}

bool Bond::isTorsionRefinable()
{
	return isNotJustForHydrogens() && !isFixed() && isUsingTorsion();
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
				vec3 pos1 = atom1->getExplicitModel()->getManyPositions()->at(0).start;
				vec3 pos2 = getMinor()->getExplicitModel()->getManyPositions()->at(0).start;
				vec3 pos3 = atom3->getExplicitModel()->getManyPositions()->at(0).start;

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
		heavy->getModel()->refreshPositions();
		_heavyAlign.lock()->getModel()->refreshPositions();
		vec3 newHPos = heavy->getModel()->getAbsolutePosition();
		vec3 origHPos = _heavyAlign.lock()->getModel()->getAbsolutePosition();
		getMinor()->getModel()->refreshPositions();
		getMajor()->getModel()->refreshPositions();
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

	addDoubleProperty("length", &_bondLength);    
	addDoubleProperty("torsion", &_torsion);    
	addDoubleProperty("exp_angle", &_expectedAngle);    
	addDoubleProperty("circle", &_circlePortion);    
	addDoubleProperty("ratio", &_geomRatio);    
	
	if (!hasWhack())
	{
		addDoubleProperty("kick", &_kick);    
	}
	else
	{
		addReference("whack", getWhack());
	}
	
	addDoubleProperty("phi", &_phi);    
	addDoubleProperty("psi", &_psi);    
	addDoubleProperty("occupancy", &_occupancy);
	addBoolProperty("occ_reset", &_resetOccupancy);
	addBoolProperty("split", &_split);

	addBoolProperty("fixed", &_fixed);
	addBoolProperty("using_torsion", &_usingTorsion);
	addBoolProperty("refine_bond_angle", &_refineBondAngle);
	addBoolProperty("refine_flexibility", &_refineFlexibility);
	addBoolProperty("disabled", &_disabled);

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

	if (category == "minor")
	{
		AtomPtr atom = ToAtomPtr(object);
		_minor = atom;
	}
	else if (category == "major")
	{
		AtomPtr atom = ToAtomPtr(object);
		_major = atom;
	}
	else if (category == "whack")
	{
		WhackPtr whack = ToWhackPtr(object);
		_whack = whack;
	}
	else if (category == "extra_sample")
	{
		AtomPtr atom = ToAtomPtr(object);
		addExtraTorsionSample(atom);
	} 
}

void Bond::postParseTidy()
{	

}

void Bond::equaliseOccupancies()
{
	double occ = 1 / (double)downstreamBondGroupCount();

	for (int i = 0; i < downstreamBondGroupCount(); i++)
	{
		for (int j = 0; j < downstreamBondCount(i); j++)
		{
			BondPtr bond = downstreamBond(i, j);
			setOccupancy(&*bond, occ);
		}
	}
	
	propagateChange();
}
