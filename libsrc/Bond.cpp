//
//  Bond.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bond.h"
#include "Atom.h"
#include "Anchor.h"
#include "fftw3d.h"
#include <iostream>
#include "Shouter.h"
#include "AtomGroup.h"
#include "Element.h"
#include "Absolute.h"
#include "maths.h"
#include <sstream>
#include <iomanip>
#include "Monomer.h"
#include "Molecule.h"
#include "Anisotropicator.h"

Bond::Bond()
{
	_anchored = false;
}

Bond::Bond(AtomPtr major, AtomPtr minor, int group)
{
	_anchored = false;
	_usingTorsion = false;
	_activated = false;
	_major = major;
	_minor = minor;
	_activeGroup = 0;
	_dampening = INITIAL_DAMPENING;
	_bendBlur = 0;
	_bondLength = 0;
	_changedPos = true;
	_changedSamples = true;
	_fixed = false;
	_blocked = false;
	_blurTotal = 0;

	BondGroup aGroup;
	aGroup.torsionAngle = 0;
	aGroup.torsionBlur = 0.0;
	aGroup.torsionVertBlur = 0.0;
	aGroup.magicAxis = make_randomish_axis();
	aGroup.magicAngle = 0;
	aGroup.occupancy = 1;
	_bondGroups.push_back(aGroup);

	_disabled = (!major || !minor);

	if (_disabled)
	{
		return;
	}

	if (getMinor()->getModel()->getClassName() == "Bond")
	{
		std::cout << "Warning!" << std::endl;
	}

	vec3 majorPos = getMajor()->getInitialPosition();
	vec3 minorPos = getMinor()->getInitialPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	_bondLength = vec3_length(difference);

	deriveBondLength();
	vec3_set_length(&difference, _bondLength);

	_bondDirection = difference;
//	_bondGroups[_activeGroup].magicAxis = _bondDirection;
	deriveBondLength();

	ModelPtr upModel = getMajor()->getModel();

	if (upModel->getClassName() == "Bond")
	{
		ToBondPtr(upModel)->addDownstreamAtom(minor, group);
	}

	if (upModel->getClassName() == "Absolute" || upModel->getClassName() == "Anchor")
	{
		ToAbsolutePtr(upModel)->addNextAtom(minor);
	}
}

Bond::Bond(Bond &other)
{
	_usingTorsion = other._usingTorsion;
	_activated = other._activated;
	_major = other._major;
	_minor = other._minor;
	_activeGroup = other._activeGroup;
	_bondGroups = other._bondGroups;
	_fixed = other._fixed;
	_disabled = other._disabled;
	_blocked = false;
	_absolute = other._absolute;
	
	_dampening = other._dampening;
	_bendBlur = 0;
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

void Bond::addDownstreamAtom(AtomPtr atom, int group, bool skipGeometry)
{
	while (_bondGroups.size() <= group)
	{
		BondGroup newGroup;
		newGroup.torsionAngle = 0;
		newGroup.torsionBlur = 0;
		newGroup.torsionVertBlur = 0.0;
		newGroup.magicAngle = 0;
		newGroup.magicAxis = make_randomish_axis();
		newGroup.occupancy = 1;
		newGroup._changedSamples = true;
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

	if (_bondGroups[group].atoms.size() > 0 && atom->getElement()->electronCount() > 1)
	{
		/* Calculate from data */
		if (!_heavyAlign.expired() && _usingTorsion)
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
	}

	AtomValue newAtom;
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

void Bond::activate(AtomGroupPtr group, AtomPtr inherit)
{
	if (_disabled)
		return;

	getMinor()->setModel(shared_from_this());

	if (group)
	{
		BondPtr myself = std::static_pointer_cast<Bond>(shared_from_this());
		group->addBond(myself);
	}

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
	mat3x3 test = mat3x3_rhbasis(xNew, reverse_bond);

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
			shout_at_helen("Torsion angle is nan!");
		}


		*newAngle = angle;
	}

	return test;
}

void Bond::setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign)
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

	makeTorsionBasis(hPos, maPos, miPos, lPos, &_bondGroups[0].torsionAngle);

	_usingTorsion = true;
}
/*
std::string Bond::getPDBContribution()
{
	AtomPtr atom = getMinor();
	std::string atomName = atom->getAtomName();
	ElementPtr element = atom->getElement();

	int tries = 10;

	if (element->getSymbol() == "H")
	{
		tries = 1;
	}

	std::ostringstream stream;
	std::vector<BondSample> *positions = getManyPositions(BondSampleThorough);

	double skip = (double)positions->size() / 25.;

	if (skip < 0) skip = 1;
	const int side = 5;
	int count = 0;

	for (double i = 0; i < positions->size(); i+= 1)
	{
		int l = i / (side * side);
		int k = (i - (l * side * side)) / side;
		int h = (i - l * side * side - k * side);

		if ((h + k) % 2 != 0 || (k + l) % 2 != 0 || (l + h) % 2 != 0)
		{
			continue;
		}

		vec3 placement = (*positions)[i].start;
		double occupancy = (*positions)[i].occupancy;

		stream << atom->pdbLineBeginning(count);
		stream << std::fixed << std::setw(8) << std::setprecision(3) << placement.x;
		stream << std::setw(8) << std::setprecision(3) << placement.y;
		stream << std::setw(8) << std::setprecision(3) << placement.z;
		stream << std::setw(6) << std::setprecision(2) << occupancy / double(tries);
		stream << std::setw(6) << std::setprecision(2) << 0;
		stream << "          ";
		stream << std::setw(2) << element->getSymbol();
		stream << "  " << std::endl;

		count++;
	}

	return stream.str();
}*/

vec3 meanOfManyPositions(std::vector<BondSample> *positions)
{
	vec3 sum = make_vec3(0, 0, 0);
	double weights = 0;

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 tmp = (*positions)[i].start;
		vec3_mult(&tmp, (*positions)[i].occupancy);
		sum = vec3_add_vec3(sum, tmp);
		weights += (*positions)[i].occupancy;
	}

	vec3_mult(&sum, 1 / weights);

	return sum;
}

std::vector<vec3> Bond::polymerCorrectedPositions()
{
	std::vector<vec3> posOnly;
	std::vector<BondSample> *positions = getManyPositions(BondSampleThorough);

	MoleculePtr molecule = getMinor()->getMolecule();
	std::vector<vec3> offsets;
	std::vector<vec3> rotationCentres;
	std::vector<mat3x3> rotations;

	if (molecule)
	{
		offsets = molecule->getCentroidOffsets();
		rotations = molecule->getRotationCorrections();
		rotationCentres = molecule->getRotationCentres();
	}

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 subtract = positions->at(i).start;

		// perform rotation element of superposition

		// remove translation aspect of superposition

	    if (rotations.size() > i && rotationCentres.size() > i)
		{
			vec3 tmp = vec3_add_vec3(subtract, rotationCentres[i]);
			mat3x3_mult_vec(rotations[i], &tmp);
			subtract = vec3_subtract_vec3(tmp, rotationCentres[i]);
		}

		if (offsets.size() > i)
		{
			subtract = vec3_subtract_vec3(subtract, offsets[i]);
		}

		posOnly.push_back(subtract);
	}

	return posOnly;
}

std::vector<BondSample> Bond::getFinalPositions()
{
	std::vector<BondSample> *positions = getManyPositions(BondSampleThorough);
	std::vector<BondSample> copyPos = *positions;
	std::vector<vec3> posOnly = polymerCorrectedPositions();

	for (int i = 0; i < copyPos.size(); i++)
	{
		if (i < posOnly.size())
		{
			copyPos[i].start = posOnly[i];
		}
	}

	_absolute = meanOfManyPositions(&copyPos);

	return copyPos;
}

FFTPtr Bond::getDistribution(bool absOnly)
{
	std::vector<BondSample> positions = getFinalPositions();

	if (absOnly) return FFTPtr();

	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2 * MAX_SCATTERING_DSTAR);

	double realLimits = (scale * n);

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);
	double occSum = 0;

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 placement = positions[i].start;

		vec3 relative = vec3_subtract_vec3(placement, _absolute);
		double occupancy = positions[i].occupancy;
		occSum += occupancy;

		vec3_mult(&relative, 1 / realLimits);

		if (relative.x != relative.x)
		{
			continue;
		}

		fft->addToReal(relative.x, relative.y, relative.z, occupancy);
	}

	fft->createFFTWplan(1);
	fft->fft(1);
	fft->invertScale();

	return std::make_shared<FFT>(*fft);
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

void Bond::calculateMagicAxis()
{
	_bondGroups[0].magicAxis = longestAxis();
	propagateChange();
}

void Bond::calculateInitialMagicAxis()
{
	BondPtr downBond = ToBondPtr(shared_from_this());
	int count = 0;
	vec3 posSum = make_vec3(0, 0, 0);
	vec3 startVec = getStaticPosition();

	while (count < FUTURE_RESIDUES && downBond->downstreamAtomGroupCount()
		   && downBond->downstreamAtomCount(0))
	{
		AtomPtr downAtom = downBond->downstreamAtom(0, 0);
		vec3 downPos = downAtom->getModel()->getStaticPosition();
		vec3 downDiff = vec3_subtract_vec3(downPos, startVec);
		vec3_set_length(&downDiff, 1);
		startVec = downPos;
		posSum = vec3_add_vec3(posSum, downDiff);

		if (downAtom->getModel()->isBond())
		{
			downBond = ToBondPtr(downAtom->getModel());
		}
		else
		{
			break;
		}

		count++;
	}

	vec3_set_length(&posSum, 1);

	if (downBond->downstreamAtomGroupCount())
	{
		_bondGroups[0].magicAxis = posSum;
	}

	propagateChange();
}

mat3x3 Bond::getMagicMat()
{
	vec3 magicAxis = _bondGroups[_activeGroup].magicAxis;
	double magicAngle = _bondGroups[_activeGroup].magicAngle;

	vec3 xAxis = make_vec3(1, 0, 0);
	vec3 zAxis = make_vec3(0, 0, 1);
	vec3 yAxis = make_vec3(0, 1, 0);

	/* Find the twizzle to put z axis onto the magic axis (around the x) */
	mat3x3 firstTwizzle = mat3x3_closest_rot_mat(magicAxis, zAxis, xAxis);

	return firstTwizzle;

	/* Find where this would place the Y axis */
	mat3x3_mult_vec(firstTwizzle, &yAxis);

	/* Find what the appropriate X axis would be */
	vec3 cross = vec3_cross_vec3(magicAxis, yAxis);

	/* Finally reconstruct the matrix from the X/Z basis vectors */
	mat3x3 magicBase = mat3x3_rhbasis(cross, magicAxis);
	mat3x3 magicRot = mat3x3_unit_vec_rotation(zAxis, magicAngle);

	mat3x3 magicInv = mat3x3_mult_mat3x3(magicRot, magicBase);

	return mat3x3_inverse(magicInv);
}

std::vector<BondSample> Bond::getCorrectedAngles(std::vector<BondSample> *prevs,
												 double circleAdd,
												 double myTorsion, double ratio)
{
	std::vector<BondSample> set;
	set.reserve(prevs->size());
	const vec3 none = make_vec3(0, 0, 0);

	_blurTotal = 0;
	AtomPtr nextAtom = downstreamAtom(_activeGroup, 0);
	vec3 myPerfectPos = getStaticPosition();
	BondPtr nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());
	vec3 nextPerfectPos = nextBond->getStaticPosition();
	double nextRatio = getGeomRatio(_activeGroup, 0);
	vec3 nextPerfectVec = vec3_subtract_vec3(nextPerfectPos, myPerfectPos);
	vec3_set_length(&nextPerfectVec, 1);

	mat3x3 magicMat = getMagicMat();

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
		vec3 nextDifference = vec3_subtract_vec3(myPerfectPos, myCurrentPos);

		/* Find out what this deviation is if beam axis is set to z */
		mat3x3_mult_vec(magicMat, &nextDifference);

		double notZ = sqrt(nextDifference.y * nextDifference.y +
						   nextDifference.x * nextDifference.x);
		double tanX = nextDifference.z / notZ;
		double xValue = sin(atan(tanX));
		double yValue = sqrt(1 - xValue * xValue);

		if (yValue != yValue)
		{
			yValue = 1;
		}

		/* We want to correct if the deviation is close to the magic angle */
		double modulation = xValue;

		if (modulation != modulation)
		{
			modulation = 0;
		}

		double rotAngle = 0;

		mat3x3_closest_rot_mat(nextBondVec, nextPerfectVec, myBondVec, &rotAngle);

		if (rotAngle != rotAngle)
		{
			rotAngle = 0;
		}

		double undoBlur = 0;
		undoBlur = rotAngle;
		undoBlur *= modulation;
		undoBlur *= fabs(_dampening);

		_blurTotal += fabs(undoBlur);

		/* This will only really apply for a kicked bond */
		double addBlur = _bondGroups[_activeGroup].torsionBlur;
		addBlur *= yValue;

		if (isFixed())
		{
			undoBlur = 0; addBlur = 0;
		}

		BondSample simple;
		simple.torsion = myTorsion + undoBlur + addBlur;
		simple.occupancy = 1;
		simple.basis = make_mat3x3();
		simple.start = nextCurrentPos;
		set.push_back(simple);
	}

	_blurTotal /= (double)prevs->size();

	return set;
}

std::vector<BondSample> *Bond::getManyPositions(BondSampleStyle style)
{
	std::vector<BondSample> *newSamples;

	if (style == BondSampleStatic)
	{
		newSamples = &_bondGroups[_activeGroup].staticSample;
	}
	else if (style == BondSampleMonteCarlo)
	{
		newSamples = &_bondGroups[_activeGroup].singleStateSample;
	}
	else
	{
		newSamples = &_bondGroups[_activeGroup].storedSamples;
	}

	bool staticAtom = (style == BondSampleStatic);

	if (!_changedSamples && style == BondSampleThorough)
	{
		return &_bondGroups[_activeGroup].storedSamples;
	}

	if (!_changedPos && style == BondSampleStatic)
	{
		return &_bondGroups[_activeGroup].staticSample;
	}

	if (_anchored && style == BondSampleStatic)
	{
		return &_bondGroups[_activeGroup].staticSample;
	}

	ModelPtr model = getMajor()->getModel();

	newSamples->clear();

	if (model->getClassName() == "Absolute")
	{
		std::vector<BondSample> *absPos = model->getManyPositions(style);
		mat3x3 magicMat = getMagicMat();

		/* We must be connected to something else, oh well */
		/* Torsion basis must be the same. */

		for (int i = 0; i < absPos->size(); i++)
		{
			vec3 majorPos = (*absPos)[i].start;
			vec3 heavyPos = getHeavyAlign()->getInitialPosition();
			vec3 none = {0, 0, 1};
			vec3 actualMajor = getMajor()->getPosition();
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

	BondPtr prevBond = std::static_pointer_cast<Bond>(model);
	prevBond = Anchor::sanitiseBond(this, prevBond);
	int myGroup = -1;
	double torsionNumber = prevBond->downstreamAtomNum(getMinor(), &myGroup);

	double totalAtoms = prevBond->downstreamAtomCount(myGroup);

	if (myGroup >= 0) // otherwise, might be next to anchor.
	{
		prevBond->setActiveGroup(myGroup);
	}

	vec3 nextPos, myMeanPos;
	double nextRatio = 0;
	BondPtr nextBond;

	if (!staticAtom && downstreamAtomCount(_activeGroup) > 0)
	{
		AtomPtr nextAtom = downstreamAtom(_activeGroup, 0);
		nextPos = nextAtom->getPosition();
		nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());
		nextRatio = getGeomRatio(_activeGroup, 0);

		myMeanPos = getStaticPosition();
	}

	std::vector<BondSample> *prevSamples = prevBond->getManyPositions(style);

	/* This is just to get a set of angles, no bases */
	double spread = 0;
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

	bool usingCompensation = (!staticAtom && isUsingTorsion()
							  && nextBond && !isFixed());

	if (usingCompensation)
	{
		myTorsions = getCorrectedAngles(prevSamples, circleAdd,
										  myTorsion, ratio);
	}
	else
	{
		myTorsions.clear();

		for (int i = 0; i < prevSamples->size(); i++)
		{
			BondSample simple;
			simple.torsion = _bondGroups[_activeGroup].torsionAngle;
			simple.basis = make_mat3x3();
			simple.start = make_vec3(0, 0, 0);
			simple.occupancy = 1;
			myTorsions.push_back(simple);
		}
	}

	for (int i = 0; i < (*prevSamples).size(); i++)
	{
		double currentTorsion = (*prevSamples)[i].torsion + circleAdd;
		double occupancy = getOccupancy(&*prevBond);

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
		(*prevSamples)[i].occupancy * occupancy;

		if (style == BondSampleStatic)
		{
			nextSample.occupancy = 1;
		}

		newSamples->push_back(nextSample);

		if (style == BondSampleStatic)
		{
			_changedPos = false;
			return newSamples;
		}
	}

	if (style == BondSampleThorough)
	{
		_changedSamples = false;
	}
	else if (style == BondSampleStatic)
	{
		_changedPos = false;
	}

	return newSamples;
}

/* This gets the position of the minor atom. */
vec3 Bond::getStaticPosition()
{
	if (!_changedPos)
	{
		return _bondGroups[_activeGroup].staticSample[0].start;
	}

	getManyPositions(BondSampleStatic);
	return _bondGroups[_activeGroup].staticSample[0].start;
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

void Bond::propagateChange()
{
	if (_anchored)
	{
		return;
	}

	_changedPos = true;
	_changedSamples = true;

	for (int j = 0; j < downstreamAtomGroupCount(); j++)
	{
		for (int i = 0; i < downstreamAtomCount(j); i++)
		{
			ModelPtr model = downstreamAtom(j, i)->getModel();
			BondPtr bond = ToBondPtr(model);
			bond->propagateChange();
		}
	}
}

ModelPtr Bond::getParentModel()
{
	AtomPtr atom = getMajor();
	ModelPtr model = atom->getModel();

	return model;
}

double Bond::getBendAngle(void *object)
{
	Bond *bond = static_cast<Bond *>(object);
	AtomPtr atom = bond->getMajor();
	ModelPtr model = atom->getModel();

	if (model->getClassName() != "Bond")
	{
		shout_at_helen("Helen should never have let this happen.\n"\
					   "Helen has tried to refine a bend connected\n"\
					   "to something that is not a bond.");
	}

	int myGroup = -1;
	BondPtr newBond = std::static_pointer_cast<Bond>(model);
	int i = newBond->downstreamAtomNum(bond->getMinor(), &myGroup);

	if (i >= 0)
	{
		double ratio = newBond->getGeomRatio(myGroup, i);
		double angle = atan(ratio) + M_PI / 2;
		return angle;
	}

	return 0;
}

void Bond::setBendAngle(void *object, double value)
{
	Bond *bond = static_cast<Bond *>(object);
	AtomPtr atom = bond->getMajor();
	ModelPtr model = atom->getModel();

	if (model->getClassName() != "Bond")
	{
		shout_at_helen("Helen should never have let this happen.\n"\
					   "Helen has tried to refine a bend connected\n"\
					   "to something that is not a bond.");
	}

	int myGroup = -1;
	BondPtr newBond = std::static_pointer_cast<Bond>(model);
	int i = newBond->downstreamAtomNum(bond->getMinor(), &myGroup);

	if (i >= 0)
	{
		double ratio = tan(value - M_PI / 2);
		newBond->setGeomRatio(myGroup, i, ratio);
	}

	newBond->propagateChange();
}

bool Bond::splitBond()
{
	BondPtr me = std::static_pointer_cast<Bond>(shared_from_this());
	BondPtr parent = std::static_pointer_cast<Bond>(getParentModel());
	int last = parent->downstreamAtomGroupCount();
	me->duplicateDownstream(parent, last);
	double occ = getOccupancy(&*parent);
	parent->setActiveGroup(last);
	double torsion = getTorsion(&*parent);

	std::cout << "Splitting bond... new torsions " << rad2deg(torsion);
	torsion += deg2rad(180);

	if (torsion > deg2rad(360))
	{
		torsion -= deg2rad(360);
	}

	std::cout << "º and " << rad2deg(torsion) << "º." << std::endl;

	setOccupancy(&*parent, occ / 2);
	setTorsion(&*parent, torsion);
	parent->setActiveGroup(0);
	setOccupancy(&*parent, occ / 2);

	propagateChange();
	return true;
}

void Bond::duplicateDownstream(BondPtr newBranch, int groupNum)
{
	ModelPtr model = getParentModel();

	if (model->getClassName() != "Bond")
	{
		return;
	}

	BondPtr myParent = std::static_pointer_cast<Bond>(model);

	BondPtr duplBond = BondPtr(new Bond(*this));
	AtomPtr duplAtom = AtomPtr(new Atom(*getMinor()));
	duplAtom->inheritParents();
	duplAtom->setModel(duplBond);
	duplBond->setMajor(newBranch->getMinor());
	duplBond->setMinor(duplAtom);

	int group = 0;
	int torsionNum = myParent->downstreamAtomNum(getMinor(), &group);
	double ratio = myParent->getGeomRatio(group, torsionNum);
	double torsion = myParent->getTorsion(group);

	newBranch->addDownstreamAtom(duplAtom, groupNum);
	int myNum = newBranch->downstreamAtomCount(groupNum) - 1;
	newBranch->setGeomRatio(groupNum, myNum, ratio);
	newBranch->setActiveGroup(groupNum);
	setTorsion(&*newBranch, torsion);
	newBranch->setActiveGroup(0);

	if (!downstreamAtomGroupCount())
	{
		return;
	}
	
	for (int i = 0; i < downstreamAtomCount(0); i++)
	{
		BondPtr nextBond = std::static_pointer_cast<Bond>(downstreamAtom(0, i)->getModel());

		if (nextBond->getClassName() != "Bond")
		{
			continue;
		}

		nextBond->duplicateDownstream(duplBond, 0);
	}
}

void Bond::setOccupancy(void *object, double value)
{
	Bond *bond = static_cast<Bond *>(object);
	bond->_bondGroups[bond->_activeGroup].occupancy = value;

	double occTotal = 0;

	for (int i = 0; i < bond->downstreamAtomGroupCount() - 1; i++)
	{
		occTotal += bond->_bondGroups[i].occupancy;
	}

	int last = bond->downstreamAtomGroupCount() - 1;
	bond->_bondGroups[last].occupancy = 1 - occTotal;

	static_cast<Bond *>(object)->propagateChange();
}

std::string Bond::shortDesc()
{
	std::ostringstream stream;
	stream << getMajor()->getMonomer()->getResidueNum() <<
	getMajor()->getAtomName() << "-" << getMinor()->getAtomName();
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
		stream << "\t" << downstreamAtom(i, 0)->getAtomName() << std::endl;
	}

	return stream.str();
}

double Bond::getMeanSquareDeviation()
{
	std::vector<BondSample> positions = getFinalPositions();

	double meanX = 0; double meanY = 0; double meanZ = 0;

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 pos = positions[i].start;
		vec3 diff = vec3_subtract_vec3(pos, _absolute);
		meanX += diff.x * diff.x;
		meanY += diff.y * diff.y;
		meanZ += diff.z * diff.z;
	}

	double score = fabs(meanX + meanY + meanZ) / 3;
	score /= (double)positions.size();
	score *= 8 * M_PI * M_PI;

	return score;
}

vec3 Bond::longestAxis()
{
	std::vector<BondSample> positions = getFinalPositions();

	std::vector<vec3> points;

	for (int i = 0; i < positions.size(); i++)
	{
		points.push_back(positions[i].start);
	}

	Anisotropicator tropicator;
	tropicator.setPoints(points);
	_realSpaceTensor = tropicator.getTensor();
	vec3 longest = tropicator.longestAxis();

	return longest;
}

mat3x3 Bond::getRealSpaceTensor()
{
	longestAxis();
	return _realSpaceTensor;
}

std::vector<AtomPtr> Bond::importantAtoms()
{
	std::vector<AtomPtr> atoms;

	for (int i = 0; i < downstreamAtomCount(_activeGroup); i++)
	{
		atoms.push_back(downstreamAtom(_activeGroup, i));
	}

	for (int i = 0; i < extraTorsionSampleCount(_activeGroup); i++)
	{
		atoms.push_back(extraTorsionSample(_activeGroup, i));
	}

	return atoms;
}

void Bond::resetBondDirection()
{
	vec3 majorPos = getMajor()->getPosition();
	vec3 minorPos = getMinor()->getPosition();

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

void Bond::setupSampling()
{
	setupGrid();
	ModelPtr model = shared_from_this();
	BondPtr bond = std::static_pointer_cast<Bond>(model);
	setJobName("bond_" + shortDesc());
	setSilent(true);

	addTorsion(bond, deg2rad(5), deg2rad(0.4));
	addSampled(importantAtoms());
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

double Bond::getFlexibilityPotential()
{
	std::vector<BondSample> *samples = getManyPositions(BondSampleThorough);
	std::vector<BondSample> *statPos = getManyPositions(BondSampleStatic);
	double sum = 0;
	double weights = 0;

	mat3x3 statBasis = statPos->at(0).basis;
	vec3 statBondDir = mat3x3_axis(statBasis, 2);

	for (int i = 0; i < samples->size(); i++)
	{
		mat3x3 iTorsionBasis = samples->at(i).basis;
		vec3 iBondDir = mat3x3_axis(iTorsionBasis, 2);

		double weight = 1;
		double angle = vec3_angle_with_vec3(iBondDir, statBondDir);

		if (angle != angle) angle = 0;

		sum += angle * weight;
		weights += weight;
	}

	double average = sum / weights;

	if (average != average)
	{
		average = 0;
	}

	double val = average;

	return val;
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

			for (int k = 0; k < downstreamAtomCount(i); k++)
			{
				AtomPtr atom3 = downstreamAtom(i, (k + 1) % downstreamAtomCount(i));

				double angle = Atom::getAngle(atom1, getMinor(), atom3);

				if (angle < 0)
				{
					/* Has no geometric entry */
					continue;
				}

				vec3 pos1 = atom1->getModel()->getStaticPosition();
				vec3 pos2 = getMinor()->getModel()->getStaticPosition();
				vec3 pos3 = atom3->getModel()->getStaticPosition();

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

				if (diff > 1e-6)
				{
					std::cout << shortDesc() << " angle "
					<< atom1->getAtomName() << "-" << getMinor()->getAtomName() << "-"
					<< atom3->getAtomName() << "\t"
					<< rad2deg(realAngle) << "\t" << rad2deg(diff) << std::endl;
				}

				ok *= (diff < 1e-6);
			}
		}
	}

	return ok;
}
