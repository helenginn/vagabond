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
	_activeGroup = 0;
	_dampening = Options::getDampen();
	_bondLength = 0;
	_changedPos = true;
	_changedSamples = true;
	_refineBondAngle = false;
	_fixed = false;
	_occupancy = 1.0;
	_occMult = 1.0;
	_torsionStepMult = 1.0;
	_anisotropyExtent = 0.0;
	_bondDirection = make_vec3(0, 0, 0);
	double initialKick = Options::getKick();

	BondGroup aGroup;
	aGroup.torsionAngle = 0;
	aGroup.torsionBlur = initialKick;
	aGroup.magicPhi = deg2rad(0);
	aGroup.magicPsi = deg2rad(0);
	_bondGroups.push_back(aGroup);
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

	ModelPtr upModel = getMajor()->getModel();

	if (upModel->getClassName() == "Bond")
	{
		ToBondPtr(upModel)->addDownstreamAtom(minor, group);
	}

	if (upModel->getClassName() == "Absolute")
	{
		ToAbsolutePtr(upModel)->addNextAtom(minor);
	}
}

Bond::Bond(Bond &other)
{
	_occupancy = other._occupancy;
	_occMult = other._occMult;
	_torsionStepMult = other._torsionStepMult;
	_usingTorsion = other._usingTorsion;
	_activated = other._activated;
	_major = other._major;
	_minor = other._minor;
	_activeGroup = other._activeGroup;
	_bondGroups = other._bondGroups;

	for (int i = 0; i < _bondGroups.size(); i++)
	{
		_bondGroups[i].atoms.clear();
	}

	_fixed = other._fixed;
	_disabled = other._disabled;
	_anchored = other._anchored;
	_absolute = other._absolute;
	_molecule = other._molecule;

	_refineBondAngle = other._refineBondAngle;
	_dampening = other._dampening;
	_bondLength = other._bondLength;
	_changedPos = true;
	_changedSamples = true;
	_bondDirection = other._bondDirection;
	_heavyAlign = other._heavyAlign;
	_lightAlign = other._lightAlign;
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
		std::cout << "Unassigned bond length for " << getMajor()->getAtomName()
		<< " to " << getMinor()->getAtomName() << "." << std::endl;
	}
}

double Bond::deriveBondAngle(AtomPtr atom)
{
	double angle = 0;
	angle = Atom::getAngle(getMajor(), getMinor(), atom);

	if (angle < 0 && atom->getElement()->electronCount() > 1)
	{
		std::cout << "Unassigned angle (" << getMajor()->getAtomName() << " to " <<
		getMinor()->getAtomName() << " to " << atom->getAtomName() << ")!" << std::endl;
	}

	return angle;
}

void Bond::resetBondAngles()
{
	for (int i = 0; i < downstreamAtomGroupCount(); i++)
	{
		for (int j = 0; j < downstreamAtomCount(i); j++)
		{
			double angle = deriveBondAngle(downstreamAtom(i, j));
			_bondGroups[i].atoms[j].expectedAngle = angle;
			double ratio = tan(angle - M_PI / 2);
			_bondGroups[i].atoms[j].geomRatio = ratio;
		}	
	}	
}

void Bond::addDownstreamAtom(AtomPtr atom, int group, bool skipGeometry)
{
	while (_bondGroups.size() <= group)
	{
		BondGroup newGroup;
		newGroup.torsionAngle = 0;
		newGroup.torsionBlur = 0;
		newGroup.magicPhi = 0;
		newGroup.magicPsi = 0;
		_bondGroups.push_back(newGroup);
	}

	vec3 pos = atom->getInitialPosition();
	vec3 start = getMinor()->getInitialPosition();
	vec3 diff = vec3_subtract_vec3(pos, start);

	/* Geometry ratio derived from model (bad) */
	double angle = vec3_angle_with_vec3(_bondDirection, diff);
	angle -= M_PI / 2;
	double ratio = tan(angle);

	double portion = -10;

	if (_bondGroups[group].atoms.size() == 0)
	{
		/* This is the first atom */
		portion = 0;
	}
	else if (_bondGroups[group].atoms.size() > 0 &&
	         atom->getElement()->electronCount() > 1)
	{
		/* Calculate from data */
		if (group == 0 && !_heavyAlign.expired() && _usingTorsion)
		{
			vec3 hPos = _heavyAlign.lock()->getInitialPosition();
			vec3 maPos = getMajor()->getInitialPosition();
			vec3 miPos = getMinor()->getInitialPosition();
			vec3 lPos = atom->getInitialPosition();

			double newAngle = 0;
			makeTorsionBasis(hPos, maPos, miPos, lPos, &newAngle);

			double oldAngle = _bondGroups[0].torsionAngle;
			double increment = newAngle - oldAngle;

			/* To be replaced with what's downstream */
			portion = increment / (2 * M_PI);
		}
		else
		{
			int atomsNow = _bondGroups[group].atoms.size();
			portion = _bondGroups[0].atoms[atomsNow].circlePortion;
		}
	}

	AtomValue newAtom;
	newAtom.placeholder = NULL;
	newAtom.atom = atom;
	newAtom.geomRatio = ratio;
	newAtom.circlePortion = portion;

	/* There is an existing atom */
	bool ok = downstreamAtomCount(group) > 0;

	angle = deriveBondAngle(atom);

	if (skipGeometry)
	{
		ok = false;
	}

	if (angle > 0)
	{
		double ratio = tan(angle - M_PI / 2);
		newAtom.geomRatio = ratio;
		newAtom.expectedAngle = angle;
	}
	else
	{
		newAtom.expectedAngle = -1;
	}

	if (ok)
	{
		/* Future atoms should be defined relative to first atom. */
		AtomPtr firstAtom = downstreamAtom(group, 0);

		/* First atom exists and is not hydrogen */
		ok *= firstAtom && (atom->getElement()->electronCount() > 1);
	}

	if (ok)
	{

		/* Geometry from minor, major, first, current atom (atom) */
		/* Bond directions/dimensions akin to unit cell! */

		AtomType central = getMinor()->getGeomType();
		AtomType preceding = getMajor()->getGeomType();
		AtomType firstDownAtom = downstreamAtom(group, 0)->getGeomType();
		AtomType newType = atom->getGeomType();

		/* Organise angles to rotate y/z around x */
		/* Which means that angle_a should match x axis */
		GeomTable table = GeomTable::getGeomTable();
		double angle_c = table.getBondAngle(preceding, central, firstDownAtom);
		double angle_b = table.getBondAngle(preceding, central, newType);
		double angle_a = table.getBondAngle(firstDownAtom, central, newType);

		if (angle_a < 0 || angle_b < 0 || angle_c < 0)
		{
			ok = false;
		}

		mat3x3 bondcell = mat3x3_from_unit_cell(1, 1, 1, rad2deg(angle_a),
		                                        rad2deg(angle_b),
		rad2deg(angle_c));

		vec3 xAxis = mat3x3_axis(bondcell, 0);
		vec3 newAtomAxis = mat3x3_axis(bondcell, 1);
		vec3 firstAtomAxis = mat3x3_axis(bondcell, 2);

		/* This angle will always come out positive */
		double increment = 0;
		mat3x3_closest_rot_mat(firstAtomAxis, newAtomAxis, xAxis, &increment);
		portion = increment / (2 * M_PI);

		/* So we must try to maintain chirality */
		double diff1 = fabs(portion - newAtom.circlePortion);
		double diff2 = fabs(-portion - newAtom.circlePortion);

		diff1 = std::min(diff1, fabs(diff1 - 1));
		diff2 = std::min(diff2, fabs(diff2 - 1));

		if (diff1 > diff2)
		{
			portion *= -1;
		}

		/* Take out unnecessary >180º circle turns */
		if ((portion - newAtom.circlePortion) > 0.5)
		{
			portion -= 1;
		}
		else if ((portion - newAtom.circlePortion) < -0.5)
		{
			portion += 1;
		}

		if (portion == portion && ok)
		{
			newAtom.circlePortion = portion;
		}
	}

	_bondGroups[group].atoms.push_back(newAtom);
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

void Bond::setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign, int groupNum)
{
	// light align can be left alone, but if neither are set, give up.
	if (_disabled || !heavyAlign || (!lightAlign && _lightAlign.expired()))
	{
		return;
	}

	_heavyAlign = heavyAlign;

	if (lightAlign)
	{
		_lightAlign = lightAlign;
	}


	/* Make torsion basis.
	* Make any starting set of angles with correct Z axis. */

	vec3 hPos = _heavyAlign.lock()->getInitialPosition();
	vec3 maPos = getMajor()->getInitialPosition();
	vec3 miPos = getMinor()->getInitialPosition();
	vec3 lPos = _lightAlign.lock()->getInitialPosition();

	makeTorsionBasis(hPos, maPos, miPos, lPos,
	                 &_bondGroups[groupNum].torsionAngle);

	_usingTorsion = true;
}


FFTPtr Bond::getDistribution(bool absOnly, int)
{
	double n = fftGridLength();
	/* Don't panic, invert scale below... this is in real space */
	double maxDStar = Options::getRuntimeOptions()->getActiveCrystalDStar();
	double scale = 1.0 / (2 * maxDStar);

	double realLimits = (scale * n);

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);
	fft->createFFTWplan(1);

	Model::addRealSpacePositions(fft, empty_vec3());

	fft->fft(1);
	fft->invertScale();

	FFTPtr newPtr;
	newPtr.reset(new FFT(*fft));
	return newPtr;
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
	mat3x3 rot = make_mat3x3();
	double phi = _bondGroups[_activeGroup].magicPhi;
	double psi = _bondGroups[_activeGroup].magicPsi;

	if (phi != 0 || psi != 0)
	{
		rot = mat3x3_rot_from_angles(phi, psi);
	}

	vec3 xAxis = make_vec3(1, 0, 0);
	vec3 zAxis = make_vec3(0, 0, 1);

	/* Find the twizzle to put z axis onto the magic axis (around the x) */
	mat3x3 firstTwizzle = mat3x3_closest_rot_mat(direction, zAxis, xAxis);
	mat3x3 multed = mat3x3_mult_mat3x3(rot, firstTwizzle);

	return multed;
}

std::vector<BondSample> Bond::getCorrectedAngles(std::vector<BondSample> *prevs,
                                                 double circleAdd,
double myTorsion, double ratio)
{
	std::vector<BondSample> set;
	set.reserve(prevs->size());
	const vec3 none = make_vec3(0, 0, 0);

	AtomPtr nextAtom = downstreamAtom(_activeGroup, 0);
	BondPtr nextBond = ToBondPtr(nextAtom->getModel());
	double nextRatio = getGeomRatio(_activeGroup, 0);

	vec3 prevMinorPosAve = make_vec3(0, 0, 0);
	vec3 myCurrentPosAve = make_vec3(0, 0, 0);
	vec3 nextCurrentPosAve = make_vec3(0, 0, 0);

	/* This loop gets average positions for the previous, current and next
	* atom positions */
	for (int i = 0; i < prevs->size(); i++)
	{
		double torsionAngle = (*prevs)[i].torsion + circleAdd;

		mat3x3 oldBasis = (*prevs)[i].basis;
		vec3 prevMinorPos = (*prevs)[i].start;
		vec3 prevHeavyPos = (*prevs)[i].old_start;

		vec3 myCurrentPos = positionFromTorsion(oldBasis, torsionAngle,
		                                        ratio, prevMinorPos);
		mat3x3 newBasis = nextBond->makeTorsionBasis(prevHeavyPos, prevMinorPos,
		                                             myCurrentPos, none);

		vec3 nextCurrentPos;
		nextCurrentPos = nextBond->positionFromTorsion(newBasis, myTorsion,
		                                               nextRatio, myCurrentPos);
		myCurrentPosAve = vec3_add_vec3(myCurrentPosAve, myCurrentPos);
		nextCurrentPosAve = vec3_add_vec3(nextCurrentPosAve, nextCurrentPos);
		prevMinorPosAve = vec3_add_vec3(prevMinorPosAve, prevMinorPos);
	}

	vec3_mult(&myCurrentPosAve, 1 / (double)prevs->size());
	vec3_mult(&prevMinorPosAve, 1 / (double)prevs->size());
	vec3_mult(&nextCurrentPosAve, 1 / (double)prevs->size());
	vec3 prevBondDir = vec3_subtract_vec3(myCurrentPosAve, prevMinorPosAve);
	vec3 averageBondDir = vec3_subtract_vec3(nextCurrentPosAve, myCurrentPosAve);
	vec3_set_length(&averageBondDir, 1.);
	vec3_set_length(&prevBondDir, 1.);
	
	vec3 crossDir = vec3_cross_vec3(averageBondDir, prevBondDir);

	mat3x3 magicMat = getMagicMat(crossDir);
	
	/* Keeps track of the average kick+dampen per bond */
	double averageModulation = 0;

	for (int i = 0; i < prevs->size(); i++)
	{
		double torsionAngle = (*prevs)[i].torsion + circleAdd;

		mat3x3 oldBasis = (*prevs)[i].basis;
		vec3 prevMinorPos = (*prevs)[i].start;
		vec3 prevHeavyPos = (*prevs)[i].old_start;

		vec3 myCurrentPos = positionFromTorsion(oldBasis, torsionAngle,
		                                        ratio, prevMinorPos);
		mat3x3 newBasis = nextBond->makeTorsionBasis(prevHeavyPos, prevMinorPos,
		                                             myCurrentPos, none);

		vec3 nextCurrentPos;
		nextCurrentPos = nextBond->positionFromTorsion(newBasis, myTorsion,
		                                               nextRatio, myCurrentPos);

		vec3 nextBondVec = vec3_subtract_vec3(nextCurrentPos, myCurrentPos);
		vec3 myBondVec = vec3_subtract_vec3(myCurrentPos, prevMinorPos);

		/* Difference between perfect and deviant position of major atom */
		vec3 nextDifference = vec3_subtract_vec3(myCurrentPos, myCurrentPosAve);

		/* Find out what this deviation is if beam axis is set to z */
		mat3x3_mult_vec(magicMat, &nextDifference);

		double notZ = sqrt(nextDifference.y * nextDifference.y +
		                   nextDifference.x * nextDifference.x);
		double tanX = nextDifference.z / notZ;
		double dampValue = sin(atan(tanX));
		double yValue = sqrt(1 - dampValue * dampValue);

		if (yValue != yValue)
		{
			yValue = 0;
		}

		yValue -= 0.5; // so the kick is applied equally in both directions.

		/* We want to correct if the deviation is close to the magic angle */
		if (dampValue != dampValue)
		{
			dampValue = 0;
		}

		double rotAngle = 0;

		/* Find the best angle for dampening */
		mat3x3_closest_rot_mat(nextBondVec, averageBondDir, myBondVec, &rotAngle);

		if (rotAngle != rotAngle)
		{
			rotAngle = 0;
		}

		double undoBlur = 0;
		undoBlur = rotAngle;
		undoBlur *= dampValue;
		undoBlur *= fabs(_dampening);

		/* This will only apply for a kicked bond */
		double addBlur = _bondGroups[_activeGroup].torsionBlur;

		addBlur *= yValue;

		if (isFixed())
		{
			undoBlur = 0; addBlur = 0;
		}
		
		double totalBlur = undoBlur + addBlur;
		
		averageModulation += totalBlur;

		BondSample simple;
		simple.torsion = myTorsion + totalBlur;
		simple.occupancy = getMultOccupancy();
		simple.basis = make_mat3x3();
		simple.start = nextCurrentPos;
		set.push_back(simple);
	}
	
	averageModulation /= (double)prevs->size();

	for (int i = 0; i < set.size(); i++)
	{
		set[i].torsion -= averageModulation;
	}

	return set;
}

std::vector<BondSample> *Bond::getManyPositions()
{
	std::vector<BondSample> *newSamples;

	newSamples = &_bondGroups[_activeGroup].storedSamples;

	if (!_changedSamples)
	{
		return &_bondGroups[_activeGroup].storedSamples;
	}

	ModelPtr model = getMajor()->getModel();

	newSamples->clear();

	if (model->getClassName() == "Absolute")
	{
		std::vector<BondSample> *absPos = model->getManyPositions();
		mat3x3 magicMat = getMagicMat(_bondDirection);

		/* We must be connected to something else, oh well */
		/* Torsion basis must be the same. */

		for (int i = 0; i < absPos->size(); i++)
		{
			vec3 majorPos = (*absPos)[i].start;
			vec3 heavyPos = getHeavyAlign()->getInitialPosition();
			vec3 none = {0, 0, 1};
			vec3 actualMajor = getMajor()->getAbsolutePosition();
			vec3 start = vec3_subtract_vec3(majorPos, _bondDirection);
			vec3 perfectStart = vec3_subtract_vec3(actualMajor, _bondDirection);

			mat3x3 newBasis = makeTorsionBasis(heavyPos, actualMajor, perfectStart, none);
			double newTorsion = _bondGroups[_activeGroup].torsionAngle;

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

			double torsionAdd = _bondGroups[_activeGroup].torsionBlur * diffValue;

			BondSample newSample;
			newSample.basis = newBasis;
			newSample.start = start;
			newSample.old_start = majorPos;
			newSample.torsion = newTorsion + torsionAdd;
			newSample.occupancy = (*absPos)[i].occupancy;
			newSamples->push_back(newSample);

		}

		return newSamples;
	}

	BondPtr prevBond = boost::static_pointer_cast<Bond>(model);
	int myGroup = -1;
	double torsionNumber = prevBond->downstreamAtomNum(getMinor(), &myGroup);

	if (myGroup < 0)
	{
		std::cout << "myGroup = " << myGroup << std::endl;
		std::cout << "Current bond: " << description() << std::endl;
		std::cout << "Previous bond: " << prevBond->description() << std::endl;
	}

	bool nextBondExists = false;
	if (_bondGroups[_activeGroup].atoms.size())
	{
		nextBondExists = true;
	}

	double totalAtoms = prevBond->downstreamAtomCount(myGroup);

	if (myGroup >= 0) // otherwise, might be next to anchor.
	{
		prevBond->setActiveGroup(myGroup);
	}

	std::vector<BondSample> *prevSamples = prevBond->getManyPositions();

	/* This is just to get a set of angles, no bases */
	double circlePortion = 0;
	double circleAdd = 0;
	if (myGroup >= 0) // otherwise, might be next to anchor.
	{
		circlePortion = prevBond->getCirclePortion(myGroup, torsionNumber);
	}

	if (circlePortion < -9) /* Likely a hydrogen */
	{
		circleAdd += deg2rad(360) * torsionNumber / totalAtoms;
	}
	else
	{
		circleAdd += deg2rad(360) * circlePortion;
	}

	double myTorsion = _bondGroups[_activeGroup].torsionAngle;
	double ratio = prevBond->getGeomRatio(myGroup, torsionNumber);

	std::vector<BondSample> myTorsions;

	bool usingCompensation = (isUsingTorsion()
	                          && nextBondExists && !isFixed());

	if (usingCompensation)
	{
		myTorsions = getCorrectedAngles(prevSamples, circleAdd,
		                                myTorsion, ratio);
	}
	else
	{
		myTorsions.clear();
		myTorsions.reserve(prevSamples->size());

		for (int i = 0; i < prevSamples->size(); i++)
		{
			BondSample simple;
			simple.torsion = _bondGroups[_activeGroup].torsionAngle;
			simple.basis = make_mat3x3();
			simple.start = make_vec3(0, 0, 0);
			simple.occupancy = getMultOccupancy();
			myTorsions.push_back(simple);
		}
	}

	double occTotal = 0;
	newSamples->reserve(prevSamples->size());

	for (int i = 0; i < (*prevSamples).size(); i++)
	{
		double currentTorsion = (*prevSamples)[i].torsion + circleAdd;
		double occupancy = getMultOccupancy();

		if (torsionNumber < 0)
		{
			shout_at_helen("Something has gone horrendously wrong\n"\
			               "in the calculation of torsion angle.\n" +
			description());
		}

		vec3 prevHeavyPos = (*prevSamples)[i].old_start;
		mat3x3 oldBasis = (*prevSamples)[i].basis;
		vec3 prevMinorPos = (*prevSamples)[i].start;

		vec3 myCurrentPos = positionFromTorsion(oldBasis, currentTorsion,
		                                        ratio, prevMinorPos);

		/* Prepping for bending */
		const vec3 none = {0, 0, 1};

		mat3x3 newBasis = makeTorsionBasis(prevHeavyPos, prevMinorPos,
		                                   myCurrentPos, none);

		BondSample nextSample;
		nextSample.basis = newBasis;
		nextSample.start = myCurrentPos;
		nextSample.old_start = prevMinorPos;
		nextSample.torsion = myTorsions[i].torsion;
		nextSample.occupancy = myTorsions[i].occupancy *
		(*prevSamples)[i].occupancy;
		occTotal += nextSample.occupancy;

		newSamples->push_back(nextSample);
	}

	if (false && getMultOccupancy() < 0.9)
	{
		std::cout << "Occ total: " << shortDesc() << " " << occTotal
		<< " instead of " << std::setprecision(4) << getMultOccupancy() << std::endl;
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

	if (downstreamAtomCount(0) == 0)
	{
		return false;
	}

	for (int i = 0; i < downstreamAtomCount(0); i++)
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

	for (int k = 0; k < propagateBonds.size(); k++)
	{
		BondPtr bond = propagateBonds[k];

		for (int j = 0; j < bond->downstreamAtomGroupCount(); j++)
		{
			for (int i = 0; i < bond->downstreamAtomCount(j); i++)
			{
				AtomPtr atom = bond->downstreamAtom(j, i);
				ModelPtr model = atom->getModel();
				BondPtr bond = ToBondPtr(model);

				propagateBonds.push_back(bond);
			}
		}

		bond->_changedPos = true;
		bond->_changedSamples = true;
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

	for (int k = 0; k < propagateBonds.size(); k++)
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
	BondPtr newBond = ToBondPtr(bond->getParentModel());

	if (!newBond || !newBond->isBond())
	{
		return;
	}

	int myGroup = -1;
	int i = newBond->downstreamAtomNum(bond->getMinor(), &myGroup);

	if (i > 0)
	{
		double angle = value / (2 * M_PI);
		newBond->_bondGroups[myGroup].atoms[i].circlePortion = angle;
	}
}

double Bond::getCirclePortion(void *object)
{
	Bond *bond = static_cast<Bond *>(object);
	BondPtr newBond = ToBondPtr(bond->getParentModel());

	if (!newBond || !newBond->isBond())
	{
		return 0;
	}

	int myGroup = -1;
	int i = newBond->downstreamAtomNum(bond->getMinor(), &myGroup);

	if (i > 0)
	{
		return newBond->_bondGroups[myGroup].atoms[i].circlePortion * 2 * M_PI;
	}

	return 0;
}

bool Bond::connectsAtom(std::string name)
{
	return (getMinor()->getAtomName() == name ||
	        getMajor()->getAtomName() == name);

}

double Bond::getBendAngle(void *object)
{
	Bond *bond = static_cast<Bond *>(object);
	BondPtr newBond = ToBondPtr(bond->getParentModel());

	if (!newBond || !newBond->isBond())
	{
		return 0;
	}

	int myGroup = -1;
	int i = newBond->downstreamAtomNum(bond->getMinor(), &myGroup);

	if (i >= 0)
	{
		double ratio = newBond->getGeomRatio(myGroup, i);
		double angle = atan(ratio) + M_PI / 2;
		return angle;
	}

	return 0;
}

double Bond::getExpectedAngle()
{
	int myGroup = -1;
	BondPtr newBond = ToBondPtr(getParentModel());
	int i = newBond->downstreamAtomNum(getMinor(), &myGroup);

	if (i >= 0)
	{
		double angle = newBond->_bondGroups[myGroup].atoms[i].expectedAngle;
		return angle;
	}

	return -1;
}

void Bond::setBendAngle(void *object, double value)
{
	Bond *bond = static_cast<Bond *>(object);
	AtomPtr atom = bond->getMajor();
	ModelPtr model = atom->getModel();

	if (model->getClassName() != "Bond")
	{
		return;
		shout_at_helen("Helen should never have let this happen.\n"\
		               "Helen has tried to refine a bend connected\n"\
		"to something that is not a bond.");
	}

	int myGroup = -1;
	BondPtr newBond = boost::static_pointer_cast<Bond>(model);
	int i = newBond->downstreamAtomNum(bond->getMinor(), &myGroup);

	if (!newBond || !newBond->isBond())
	{
		return;
	}

	if (i >= 0)
	{
		double ratio = tan(value - M_PI / 2);
		newBond->setGeomRatio(myGroup, i, ratio);
	}

	newBond->propagateChange(10);
}

bool Bond::splitBond(int start)
{
	BondPtr me = boost::static_pointer_cast<Bond>(shared_from_this());
	BondPtr parent = boost::static_pointer_cast<Bond>(getParentModel());
	int last = parent->downstreamAtomGroupCount();
	BondPtr dupl = me->duplicateDownstream(parent, last, start);
	double torsion = getTorsion(&*me);
	setTorsion(&*me, torsion);

	_occupancy /= 2;
	dupl->_occupancy /= 2;


	propagateChange();

	return true;
}

BondPtr Bond::duplicateDownstream(BondPtr newBranch, int groupNum, int start)
{
	ModelPtr model = getParentModel();

	if (!isBond())
	{
		return BondPtr();
	}

	BondPtr myParent = boost::static_pointer_cast<Bond>(model);

	BondPtr duplBond = BondPtr(new Bond(*this));

	AtomPtr duplAtom = AtomPtr();
	AtomList list = getMinor()->getMonomer()->findAtoms(getMinor()->getAtomName());

	/* Do we have something in the list which has an Absolute model?
	* if not, leave as default (new atom) */

	bool changed = false;

	/* Existing atoms: list.size() */
	for (int i = 0; i < list.size(); i++)
	{
		if (list[i].expired()) continue;
		AtomPtr atom = list[i].lock();
		ModelPtr model = atom->getModel();

		if (model->isAbsolute())
		{
			duplAtom = atom;
			changed = true;
			break;
		}
	}

	if (!changed)
	{
		duplAtom = AtomPtr(new Atom(*getMinor()));
		char conformer[] = "a";
		getMinor()->setAlternativeConformer(conformer);

		for (int i = 0; i < list.size(); i++)
		{
			conformer[0]++;
		}

		duplAtom->setAlternativeConformer(conformer);
		duplAtom->setFromPDB(false);
	}

	duplAtom->inheritParents();
	duplAtom->setModel(duplBond);
	duplBond->setMajor(newBranch->getMinor());
	duplBond->setMinor(duplAtom);

	int group = 0;

	double torsion = myParent->getTorsion(group);

	/* Need to add the downstream atoms from group 0 which are not duplAtom */

	if (groupNum != 0)
	{
		for (int i = 0; i < newBranch->downstreamAtomCount(0); i++)
		{
			if (newBranch->downstreamAtom(0, i)->getAtomName() != duplAtom->getAtomName())
			{
				newBranch->addDownstreamAtom(newBranch->downstreamAtom(0, i), groupNum);
			}
		}
	}

	newBranch->addDownstreamAtom(duplAtom, groupNum);
	newBranch->setActiveGroup(groupNum);
	setTorsion(&*newBranch, torsion);
	newBranch->setActiveGroup(0);

	for (int i = 0; i < newBranch->downstreamAtomCount(groupNum); i++)
	{
		double portion = myParent->_bondGroups[0].atoms[i].circlePortion;
		newBranch->_bondGroups[groupNum].atoms[i].circlePortion = portion;
	}

	if (!downstreamAtomGroupCount())
	{
		return duplBond;
	}

	for (int i = start; i < downstreamAtomCount(0); i++)
	{
		BondPtr nextBond = boost::static_pointer_cast<Bond>(downstreamAtom(0, i)->getModel());

		if (!nextBond->isBond())
		{
			continue;
		}

		nextBond->duplicateDownstream(duplBond, 0, 0);
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
	<< rad2deg(_bondGroups[0].torsionAngle) << std::endl;
	stream << "Bond downstream groups: ("
	<< downstreamAtomGroupCount() << "):" << std::endl;
	stream << "Bond downstream atoms (first) ("
	<< downstreamAtomCount(0) << "):" << std::endl;

	for (int i = 0; i < downstreamAtomGroupCount(); i++)
	{
		for (int j = 0; j < downstreamAtomCount(i); j++)
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

void Bond::reverseDownstreamAtoms(int group)
{
	std::vector<AtomValue> newAtoms;
	newAtoms.push_back(_bondGroups[group].atoms[0]);

	for (int i = downstreamAtomCount(group) - 1; i > 0; i--)
	{
		newAtoms.push_back(_bondGroups[group].atoms[i]);
	}

	_bondGroups[group].atoms = newAtoms;
}

ModelPtr Bond::reverse(BondPtr upstreamBond)
{
	ModelPtr nextBond = getMajor()->getModel();
	AtomPtr minor = getMinor();
	AtomPtr major = getMajor();
	setMinor(major);
	setMajor(minor);

	AtomPtr heavy = getHeavyAlign();
	AtomPtr light = getLightAlign();
	_heavyAlign = light;
	_lightAlign = heavy;

	vec3_mult(&_bondDirection, -1);

	for (int i = 0; i < downstreamAtomGroupCount(); i++)
	{
		if (upstreamBond)
		{
			_bondGroups[i].atoms.erase(_bondGroups[i].atoms.begin());
			upstreamBond->_bondGroups[i].atoms.clear();

			// FIXME
			upstreamBond->addDownstreamAtom(major, i);

			for (int j = downstreamAtomCount(i) - 1; j >= 0; j--)
			{
				AtomPtr atom = downstreamAtom(i, j);
				upstreamBond->addDownstreamAtom(atom, i);
			}
		}
	}

	setActiveGroup(0);
	activate();

	return nextBond;
}

bool Bond::isRefinable()
{
	return isNotJustForHydrogens() && !isFixed() && isUsingTorsion();
}

bool Bond::test()
{
	bool ok = true;

	/* Test of geometry for multiple downstream atoms */
	for (int i = 0; i < downstreamAtomGroupCount(); i++)
	{
		for (int j = -1; j < downstreamAtomCount(i); j++)
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

			for (int k = 0; k < downstreamAtomCount(i); k++)
			{
				AtomPtr atom3 = downstreamAtom(i, (k + 1) % downstreamAtomCount(i));

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

	for (int i = 0; i < samples->size(); i++)
	{
		total += samples->at(i).occupancy;
	}

	return total;
}

bool parseAtomValue(char **blockPtr, AtomValue *atomv)
{
	//    std::cout << "Parsing atom value" << std::endl;
	char *block = *blockPtr;
	char *keyword = NULL;
	char *value = NULL;

	if (block[0] == '}')
	{
		*blockPtr = block + 1;
		incrementIndent(blockPtr);
		return false;
	}

	*blockPtr = keywordValue(block, &keyword, &value);

	double val = strtod(value, NULL);

	if (strcmp(keyword, "geom_ratio") == 0)
	{
		atomv->geomRatio = val;
	}
	else if (strcmp(keyword, "expected") == 0)
	{
		atomv->expectedAngle = val;
	} 
	else if (strcmp(keyword, "circle_add") == 0)
	{
		atomv->circlePortion = val;
	}
	else if (strcmp(keyword, "atom") == 0)
	{
		atomv->placeholder = new std::string(value);
	}

	//    printf("%s ... %s\n", keyword, value);

	block++;
	incrementIndent(&block);

	return true;
}


bool parseBondGroup(char **blockPtr, BondGroup *group)
{
	//    std::cout << "Parsing bond group" << std::endl;
	char *block = *blockPtr;
	char *keyword = NULL;
	char *value = NULL;

	block = keywordValue(block, &keyword, &value);

	if (block == NULL)
	{
		std::cout << "Bad things happened" << std::endl;
		return false;
	}

	if (strcmp(keyword, "object") == 0)
	{
		// do a thing to load an atom value.
		AtomValue atomv;

		if (block[0] != '{')
		{
			std::cout << "Was expecting a {" << std::endl;
			return false;
		}

		block++;
		incrementIndent(&block);

		bool another = true;

		while (another)
		{
			another = parseAtomValue(&block, &atomv);
		}

		if (block == NULL)
		{
			return false;
		}        

		group->atoms.push_back(atomv);
	}
	else
	{
		double val = strtod(value, NULL);

		if (strcmp(keyword, "torsion") == 0)
		{
			group->torsionAngle = val;
		}
		else if (strcmp(keyword, "kick") == 0)
		{
			group->torsionBlur = val;
		}
		else if (strcmp(keyword, "phi") == 0)
		{
			group->magicPhi = val;
		}
		else if (strcmp(keyword, "psi") == 0)
		{
			group->magicPsi = val;
		}

		//        printf("%s ... %s\n", keyword, value);
	}

	*blockPtr = block;

	if (block[0] == '}')
	{
		*blockPtr = block + 1;
		incrementIndent(blockPtr);
		return false;
	}

	return true;

}

char *Bond::decodeBondGroup(void *bond, void *bondGroup, char *block)
{
	char *start = block;
	std::vector<BondGroup> *bondGroups = NULL;
	bondGroups = static_cast<std::vector<BondGroup> *>(bondGroup); 
	
	bondGroups->clear();

	// poised at "object", hopefully. Let's check.
	while (true)
	{    
		char *white = strchrwhite(block);
		
		*white = 0;

		if (block[0] == '}')
		{
			break;
		}

		if (strcmp(block, "object") != 0)
		{
			std::cout << "Was expecting an object for a BondGroup." << std::endl;
			
			std::cout << "Instead, " << block[0] << std::endl;

			return NULL;
		}

		block = white + 1;
		incrementIndent(&block);

		/* Starting contents of object */

		if (block[0] != '{')
		{
			std::cout << "Was expecting a {" << std::endl;
			return NULL;
		}

		block++;
		incrementIndent(&block);

		bool another = true;

		BondGroup group;

		while (another)
		{
			another = parseBondGroup(&block, &group);
		}

		bondGroups->push_back(group);

		if (block == NULL)
		{
			return NULL;
		}

	}

	if (block[0] == 0) return NULL;

	//    std::cout << "Leaving on character " << block[0] << std::endl;

	return block;
}

void Bond::encodeBondGroup(void *bond, void *bondGroup,
                           std::ostream &stream, int in)
{
	std::vector<BondGroup> *groups = NULL;
	groups = static_cast<std::vector<BondGroup> *>(bondGroup);

	stream << std::setprecision(8);

	for (int i = 0; i < groups->size(); i++)
	{
		stream << indent(in) << "object " << std::endl;
		stream << indent(in) << "{" << std::endl;
		in++;
		stream << indent(in) << "torsion = " << (*groups)[i].torsionAngle << std::endl;
		stream << indent(in) << "kick = " << (*groups)[i].torsionBlur << std::endl;
		stream << indent(in) << "phi = " << (*groups)[i].magicPhi << std::endl;
		stream << indent(in) << "psi = " << (*groups)[i].magicPsi << std::endl;

		for (int j = 0; j < (*groups)[i].atoms.size(); j++)
		{
			AtomValue *atom = &(*groups)[i].atoms[j];
			stream << indent(in) << "object " << std::endl;
			stream << indent(in) << "{" << std::endl;
			in++;
			stream << indent(in) << "geom_ratio = " << atom->geomRatio << std::endl;
			stream << indent(in) << "expected = " << atom->expectedAngle << std::endl;
			stream << indent(in) << "circle_add = " << atom->circlePortion << std::endl;
			AtomPtr refAtom = atom->atom.lock();
			stream << indent(in) << "atom = " << refAtom->getAbsolutePath() << std::endl;
			in--;
			stream << indent(in) << "}" << std::endl;
		}

		in--;
		stream << indent(in) << "}" << std::endl;
	}
}

void Bond::addProperties()
{
	addReference("minor", getMinor());
	addReference("major", getMajor());
	addReference("heavy", getHeavyAlign());
	addReference("light", getLightAlign());

	addDoubleProperty("length", &_bondLength);    
	addDoubleProperty("dampening", &_dampening);
	addDoubleProperty("occupancy", &_occupancy);
	addDoubleProperty("occ_mult", &_occMult);

	addBoolProperty("fixed", &_fixed);
	addBoolProperty("anchored", &_anchored);
	addBoolProperty("using_torsion", &_usingTorsion);
	addBoolProperty("refine_bond_angle", &_refineBondAngle);
	addBoolProperty("activated", &_activated);
	addBoolProperty("disabled", &_disabled);

	addVec3Property("bond_direction", &_bondDirection);

	addCustomProperty("bond_group", &_bondGroups, this,
	                  encodeBondGroup, decodeBondGroup);

	Model::addProperties();
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
}

void Bond::postParseTidy()
{	
	/* Get real references to the downstream atoms */
	for (int i = 0; i < downstreamAtomGroupCount(); i++)
	{
		for (int j = 0; j < downstreamAtomCount(i); j++)
		{
			if (!_bondGroups[i].atoms[j].placeholder)
			{
				continue;
			}

			std::string str = *_bondGroups[i].atoms[j].placeholder;

			ParserPtr parser = Parser::resolveReference(str);
			AtomPtr atom = ToAtomPtr(parser);
			
			if (!atom)
			{
				shout_at_helen("Post parse tidy fail in Bond.");
			}

			_bondGroups[i].atoms[j].atom = atom;
			delete _bondGroups[i].atoms[j].placeholder;
			_bondGroups[i].atoms[j].placeholder = NULL;
		}
	}
	
	/* If fixed, add oneself to the last non-fixed bond as
	* 	an extra torsion sample */
	/* This may be buggy. */
	
	if (!isFixed())
	{
		return;	
	}
	
	ModelPtr lastModel = getParentModel();
	if (!lastModel) return;

	BondPtr lastBond = ToBondPtr(lastModel);
	
	while (lastBond->isFixed())
	{
		lastModel = lastBond->getParentModel();
		
		if (!lastModel) return;
		
		lastBond = ToBondPtr(lastModel);
		
		if (!lastBond) return;
	}
	
	if (lastBond->downstreamAtomGroupCount() == 0)
	{
		return;	
	}

	/* FIXME: may not be zero */
	lastBond->addExtraTorsionSample(getMinor(), 0);
}
