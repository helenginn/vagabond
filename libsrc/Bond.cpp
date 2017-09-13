//
//  Bond.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bond.h"
#include "Atom.h"
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

mat3x3 Bond::magicAxisChecks[] =
{
	{ 1, 0, 0, 0, 1, 0, 0, 0, 1},
	{-1, 0, 0, 0, 1, 0, 0, 0, 1},
	{0, -1, 0, 1, 0, 0, 0, 0, 1},
	{0, 1, 0, -1, 0, 0, 0, 0, 1},
	{0, 0, 1, 0, 1, 0, -1, 0, 0},
	{0, 0, -1, 0, 1, 0, 1, 0, 0},
};

Bond::Bond(AtomPtr major, AtomPtr minor, int group)
{
	_usingTorsion = false;
	_activated = false;
	_major = major;
	_minor = minor;
	_activeGroup = 0;
	_dampening = 0.05;
	_bendBlur = 0;
	_bondLength = 0;
	_changedPos = true;
	_changedSamples = true;
	_fixed = false;
	_blocked = false;
	_currentCheck = 0;

	BondGroup aGroup;
	aGroup.torsionAngle = 0;
	aGroup.torsionBlur = 0.0;
	aGroup.torsionVertBlur = 0.0;
	aGroup.magicAxis = make_randomish_axis();
	aGroup.hRot = 0;
	aGroup.kRot = 0;
	aGroup.occupancy = 1;
	_bondGroups.push_back(aGroup);

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
	_bondGroups[_activeGroup].magicAxis = _bondDirection;
	deriveBondLength();

	ModelPtr upModel = getMajor()->getModel();

	if (upModel->getClassName() == "Bond")
	{
		std::static_pointer_cast<Bond>(upModel)->addDownstreamAtom(minor, group);
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
	_blocked = false;
	_currentCheck = 0;

	_dampening = -0;
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

void Bond::deriveBondAngle(int group, int n)
{
	AtomType type1 = getMajor()->getGeomType();
	AtomType type2 = getMinor()->getGeomType();

	if (downstreamAtomGroupCount() <= group)
	{
		return;
	}

	AtomPtr next = downstreamAtom(group, n);
	AtomType type3 = next->getGeomType();

	GeomTable table = GeomTable::getGeomTable();
	double angle = table.getBondAngle(type1, type2, type3);

	if (angle < 0 && next->getElement()->electronCount() > 1)
	{
		std::cout << "Unassigned angle (" << getMajor()->getAtomName() << " to " <<
		getMinor()->getAtomName() << " to " << next->getAtomName() << ")!" << std::endl;
	}

	if (angle > 0)
	{
		double ratio = tan(angle - M_PI / 2);
		setGeomRatio(group, n, ratio);
	}
}

void Bond::addDownstreamAtom(AtomPtr atom, int group)
{
	while (_bondGroups.size() <= group)
	{
		BondGroup newGroup;
		newGroup.torsionAngle = 0;
		newGroup.torsionBlur = 0;
		newGroup.torsionVertBlur = 0.0;
		newGroup.hRot = 0;
		newGroup.kRot = 0;
		newGroup.magicAxis = make_randomish_axis();
		newGroup.occupancy = 1;
		newGroup._changedSamples = true;
		_bondGroups.push_back(newGroup);
	}

	vec3 pos = atom->getInitialPosition();
	vec3 start = getMinor()->getInitialPosition();
	vec3 diff = vec3_subtract_vec3(pos, start);

	double angle = vec3_angle_with_vec3(_bondDirection, diff);
	angle -= M_PI / 2;
	double ratio = tan(angle);
	double portion = -10;

	if (_bondGroups[group].atoms.size() > 0)
	{
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

			

			portion = increment / (2 * M_PI);
		}
	}

	AtomValue newAtom;
	newAtom.atom = atom;
	newAtom.geomRatio = ratio;
	newAtom.circlePortion = portion;
	_bondGroups[group].atoms.push_back(newAtom);

	int totalAtoms = _bondGroups[group].atoms.size();
	int totalGroups = downstreamAtomGroupCount();

	deriveBondAngle(totalGroups - 1, totalAtoms - 1);
}

void Bond::setMinor(AtomPtr newMinor)
{
	_minor = newMinor;

	if (_bondLength < 0)
	{
		deriveBondLength();
	}

	vec3 majorPos = getMajor()->getPosition();
	vec3 minorPos = getMinor()->getInitialPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	vec3_set_length(&difference, _bondLength);

	_bondDirection = difference;
}

void Bond::activate(AtomGroupPtr group, AtomPtr inherit)
{
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

		*newAngle = angle;
	}

	_torsionBasis = test;
	return test;
}

void Bond::setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign)
{
	_heavyAlign = heavyAlign;
	_lightAlign = lightAlign;

	/* Make torsion basis.
	 * Make any starting set of angles with correct Z axis. */

	vec3 hPos = _heavyAlign.lock()->getInitialPosition();
	vec3 maPos = getMajor()->getInitialPosition();
	vec3 miPos = getMinor()->getInitialPosition();
	vec3 lPos = _lightAlign.lock()->getInitialPosition();

	makeTorsionBasis(hPos, maPos, miPos, lPos, &_bondGroups[0].torsionAngle);

	_usingTorsion = true;
}

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
}

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

FFTPtr Bond::getDistribution()
{
	std::vector<BondSample> *positions = getManyPositions(BondSampleThorough);
	_absolute = meanOfManyPositions(positions);

	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2 * MAX_SCATTERING_DSTAR);

	double realLimits = (scale * n);

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);
	double occSum = 0;

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 placement = (*positions)[i].start;
		vec3 relative = vec3_subtract_vec3(placement, _absolute);
		double occupancy = (*positions)[i].occupancy;
		occSum += occupancy;

		vec3_mult(&relative, 1 / realLimits);

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

std::vector<BondSample> Bond::sampleMyAngles(double angle, double sigma,
											 bool singleState)
{
	double sum = 0;
	sigma = fabs(sigma);
	std::vector<BondSample> samples;

	if (singleState)
	{
		BondSample sample;
		sample.occupancy = 1;
		sample.torsion = random_norm_dist(angle, sigma);
		samples.push_back(sample);
		return samples;
	}

	if (sigma <= 0)
	{
		BondSample sample;
		sample.occupancy = 1;
		sample.torsion = angle;
		samples.push_back(sample);

		return samples;
	}

	double interval = ANGLE_SAMPLING;

	if (sigma <= interval)
	{
		interval = sigma * 1.0;
	}

	samples.reserve(sigma * 2 / interval + 1);
	int count = 0;

	for (double ang = -2.0 * sigma; ang <= 2.0 * sigma; ang += interval)
	{
		double relFreq = normal_distribution(ang, sigma);
		sum += relFreq;
		BondSample sample;
		sample.occupancy = relFreq;
		sample.torsion = ang + angle;
		samples.push_back(sample);
		count++;
	}

	for (int i = 0; i < samples.size(); i++)
	{
		samples[i].occupancy /= sum;
	}

	return samples;
}

vec3 Bond::getFixedAxis(vec3 axis, double hRot, double kRot)
{
	vec3 xAxis = make_vec3(1, 0, 0);
	vec3 yAxis = make_vec3(0, 1, 0);
	mat3x3 smallRotH = mat3x3_unit_vec_rotation(xAxis, hRot);
	mat3x3 smallRotK = mat3x3_unit_vec_rotation(yAxis, kRot);
	mat3x3 smallRots = mat3x3_mult_mat3x3(smallRotK, smallRotH);
	mat3x3 magicMat = mat3x3_mult_mat3x3(smallRots, magicAxisChecks[_currentCheck]);

	vec3 compAxis = axis;
	mat3x3_mult_vec(magicMat, &compAxis);
	_currentCheck = 0;

	return compAxis;
}

void Bond::resetAxis()
{
	vec3 compAxis = getFixedAxis(_bondGroups[_activeGroup].magicAxis,
								 _bondGroups[_activeGroup].hRot,
								 _bondGroups[_activeGroup].kRot);
	_bondGroups[_activeGroup].magicAxis = compAxis;
	_bondGroups[_activeGroup].hRot = 0;
	_bondGroups[_activeGroup].kRot = 0;
}

mat3x3 Bond::getMagicMat()
{
	vec3 xAxis = make_vec3(1, 0, 0);
	vec3 zAxis = make_vec3(0, 0, 1);

	vec3 compAxis = getFixedAxis(_bondGroups[_activeGroup].magicAxis,
								 _bondGroups[_activeGroup].hRot,
								 _bondGroups[_activeGroup].kRot);
	double compLength = vec3_length(compAxis);
	mat3x3 compMat = mat3x3_closest_rot_mat(xAxis, compAxis, zAxis);
	mat3x3_scale(&compMat, compLength, compLength, compLength);
	mat3x3 invComp = mat3x3_inverse(compMat);
	mat3x3 magicMat = mat3x3_mult_mat3x3(invComp, magicAxisChecks[_currentCheck]);

	return magicMat;
}

std::vector<BondSample> Bond::getCorrectedAngles(std::vector<BondSample> *prevs,
												 double circleAdd,
												 double myTorsion, double ratio,
												 double *transferBlur)
{
	std::vector<BondSample> set;
	BondPtr prevBond = std::static_pointer_cast<Bond>(getParentModel());
	const vec3 none = make_vec3(0, 0, 0);
	*transferBlur = 0;

	AtomPtr nextAtom = downstreamAtom(_activeGroup, 0);
	vec3 myPerfectPos = getStaticPosition();
	BondPtr nextBond = std::static_pointer_cast<Bond>(nextAtom->getModel());
	vec3 nextPerfectPos = nextBond->getStaticPosition();
	double nextRatio = getGeomRatio(_activeGroup, 0);

	const bool useMagicAxis = true;

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
		vec3 nextPerfectVec = vec3_subtract_vec3(nextPerfectPos, myPerfectPos);
		vec3 nextDifference = vec3_subtract_vec3(myPerfectPos, myCurrentPos);
		mat3x3_mult_vec(magicMat, &nextDifference);
		double notX = sqrt(nextDifference.y * nextDifference.y +
						   nextDifference.z * nextDifference.z);
		double notZ = sqrt(nextDifference.y * nextDifference.y +
						   nextDifference.x * nextDifference.x);
		double tanX = nextDifference.x / notX;
		double tanZ = nextDifference.z / notZ;

		double xValue = sin(atan(tanX));
		double yValue = cos(atan(tanX));
		double zValue = sin(atan(tanZ));

		if (xValue != xValue)
		{
			xValue = 0;
		}

		if (!useMagicAxis)
		{
			xValue = 1;
		}

		if (yValue != yValue)
		{
			yValue = 1;
		}

		if (zValue != zValue)
		{
			zValue = 0;
		}

		vec3_set_length(&myBondVec, 1);
		vec3_set_length(&nextPerfectVec, 1);

		double rotAngle = 0;
		mat3x3_closest_rot_mat(nextBondVec, nextPerfectVec, myBondVec, &rotAngle);

		if (rotAngle != rotAngle)
		{
			rotAngle = 0;
		}

		double undoBlur = 0;
		undoBlur = rotAngle;
		undoBlur *= xValue;
		undoBlur *= _dampening;

		if (_dampening <= 0) undoBlur = 0;

		double addBlur = _bondGroups[_activeGroup].torsionBlur;
		addBlur *= yValue;

		double bigBlur = addBlur;
		bigBlur *= -_dampening;
		if (_dampening > 0) bigBlur = 0;

		BondSample simple;
		simple.torsion = myTorsion + undoBlur + bigBlur;
		simple.occupancy = 1;
		simple.basis = make_mat3x3();
		simple.start = nextCurrentPos;
		set.push_back(simple);
	}

	*transferBlur /= set.size();
	*transferBlur = sqrt(*transferBlur);

	return set;
}

std::vector<BondSample> *Bond::getManyPositions(BondSampleStyle style)
{
	ModelPtr model = getMajor()->getModel();

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
	bool monteCarlo = (style == BondSampleMonteCarlo);

	if (!_changedSamples && style == BondSampleThorough)
	{
		return &_bondGroups[_activeGroup].storedSamples;
	}

	newSamples->clear();

	if (model->getClassName() == "Absolute")
	{
		std::vector<BondSample> *absPos;
		absPos = ToAbsolutePtr(model)->getManyPositions(style);
		mat3x3 magicMat = getMagicMat();

		/* We must be connected to something else, oh well */
		/* Torsion basis must be the same. */

		for (int i = 0; i < absPos->size(); i++)
		{
			vec3 majorPos = (*absPos)[i].start;
			vec3 heavyPos = getHeavyAlign()->getInitialPosition();
			vec3 none = {0, 0, 1};
			vec3 actualMajor = getMajor()->getPosition();
			vec3 halfDir = _bondDirection;
			vec3_mult(&halfDir, 0.5);
			vec3 midPoint = vec3_subtract_vec3(actualMajor, halfDir);

			vec3 newDir = vec3_subtract_vec3(midPoint, majorPos);
			vec3_mult(&newDir, 2.0);

			vec3 start = vec3_add_vec3(majorPos, newDir);
			start = vec3_subtract_vec3(majorPos, _bondDirection);


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

	std::vector<BondSample> myBendings;

	/* This is just to get a set of angles, no bases */
	double spread = staticAtom ? 0 : _bendBlur;
	myBendings = sampleMyAngles(0, spread, monteCarlo);

	double circleAdd = 0;
	double circlePortion = prevBond->getCirclePortion(myGroup, torsionNumber);

	if (circlePortion > -9)
	{
		circleAdd += deg2rad(360) * circlePortion;
	}
	else
	{
		circleAdd += deg2rad(360) * torsionNumber / totalAtoms;
	}

	double myTorsion = _bondGroups[_activeGroup].torsionAngle;
	double ratio = prevBond->getGeomRatio(myGroup, torsionNumber);

	std::vector<BondSample> myTorsions;
	std::vector<BondSample> myDupPrevSamples;
	std::vector<BondSample> tempTorsions;

	bool usingCompensation = (!staticAtom && isUsingTorsion()
							  && nextBond && !isFixed());

	if (usingCompensation)
	{
		double transferBlur = 0;
		tempTorsions = getCorrectedAngles(prevSamples, circleAdd, myTorsion,
										  ratio, &transferBlur);

		for (int i = 0; i < tempTorsions.size(); i++)
		{
			std::vector<BondSample> moreTorsions;
			moreTorsions = sampleMyAngles(tempTorsions[i].torsion, 0,
										  monteCarlo);

			for (int j = 0; j < moreTorsions.size(); j++)
			{
				myDupPrevSamples.push_back((*prevSamples)[i]);
			}

			myTorsions.reserve(myTorsions.size() + moreTorsions.size());
			myTorsions.insert(myTorsions.end(), moreTorsions.begin(), moreTorsions.end());
		}
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
			myDupPrevSamples.push_back((*prevSamples)[i]);
		}
	}

	for (int i = 0; i < myDupPrevSamples.size(); i++)
	{
		double currentTorsion = myDupPrevSamples[i].torsion + circleAdd;

		spread = staticAtom ? 0 : _bondGroups[_activeGroup].torsionBlur;

		double occupancy = getOccupancy(&*prevBond);

		if (torsionNumber < 0)
		{
			shout_at_helen("Something has gone horrendously wrong\n"\
						   "in the calculation of torsion angle.\n" +
						   description());
		}

		vec3 prevHeavyPos = myDupPrevSamples[i].old_start;
		mat3x3 oldBasis = myDupPrevSamples[i].basis;
		vec3 prevMinorPos = myDupPrevSamples[i].start;
		vec3 myCurrentPos = positionFromTorsion(oldBasis, currentTorsion,
												ratio, prevMinorPos);

		/* Prepping for bending */
		vec3 zAxis = {0, 0, 1};
		mat3x3_mult_vec(oldBasis, &zAxis);
		const vec3 none = {0, 0, 1};

		mat3x3 newBasis = makeTorsionBasis(prevHeavyPos, prevMinorPos,
										   myCurrentPos, none);

		/* do the bendy thing */

		for (int j = 0; j < myBendings.size(); j++)
		{
			vec3 cross = vec3_cross_vec3(myCurrentPos, zAxis);
			vec3_set_length(&cross, 1);
			mat3x3 bend = mat3x3_unit_vec_rotation(cross, myBendings[j].torsion);
			mat3x3_mult_vec(bend, &myCurrentPos);

			BondSample nextSample;
			nextSample.basis = newBasis;
			nextSample.start = myCurrentPos;
			nextSample.old_start = prevMinorPos;
			nextSample.torsion = myTorsions[i].torsion;
			nextSample.occupancy = myTorsions[i].occupancy *
			myBendings[j].occupancy * myDupPrevSamples[i].occupancy * occupancy;

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
	_changedPos = true;
	_changedSamples = true;

	for (int j = 0; j < downstreamAtomGroupCount(); j++)
	{
		for (int i = 0; i < downstreamAtomCount(j); i++)
		{
			ModelPtr model = downstreamAtom(j, i)->getModel();
			BondPtr bond = std::static_pointer_cast<Bond>(model);
			bond->propagateChange();
		}
	}
}

std::string Bond::description()
{
	std::ostringstream stream;
	stream << "Bond: connects " << getMajor()->getAtomName() << " to "
	<< getMinor()->getAtomName() << std::endl;
	stream << "Bond length: " << _bondLength << " Å" << std::endl;
	stream << "Bond torsion angle: " << rad2deg(_bondGroups[0].torsionAngle) << std::endl;
	stream << "Bond downstream groups: (" << downstreamAtomGroupCount() << "):" << std::endl;
	stream << "Bond downstream atoms (first) (" << downstreamAtomCount(0) << "):" << std::endl;

	for (int i = 0; i < downstreamAtomCount(0); i++)
	{
		stream << "\t" << downstreamAtom(i, 0)->getAtomName() << std::endl;
	}

	return stream.str();
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
	stream << getMajor()->getAtomName() << "-" << getMinor()->getAtomName();
	return stream.str();
}

double Bond::getMeanSquareDeviation(double target, int index)
{
	double targMult = 3;
	target /= 8 * M_PI * M_PI;
    target *= targMult;

	std::vector<BondSample> *positions = getManyPositions(BondSampleThorough);
	vec3 mean = make_vec3(0, 0, 0);

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 pos = (*positions)[i].start;

		mean = vec3_add_vec3(mean, pos);
	}

	double mult = 1 / (double)positions->size();
	vec3_mult(&mean, mult);

	double meanSq = 0;
	double meanX = 0; double meanY = 0; double meanZ = 0;

	for (int i = 0; i < positions->size(); i++)
	{
		vec3 pos = (*positions)[i].start;
		double occupancy = sqrt((*positions)[i].occupancy);

		vec3 diff = vec3_subtract_vec3(pos, mean);
		vec3_mult(&diff, occupancy);
		meanX += diff.x * diff.x;
		meanY += diff.y * diff.y;
		meanZ += diff.z * diff.z;
		meanSq += vec3_sqlength(diff);
	}

	double score = 0;

	if (target >= 0)
	{
		score += fabs(meanX - target / 3);
		score += fabs(meanY - target / 3);
		score += fabs(meanZ - target / 3);
	}
	else
	{
		score = fabs(meanSq - target);
		score *= 8 * M_PI * M_PI;
		score /= targMult;
	}

	if (index >= 0)
	{
		double val = 0;
		switch (index) {
			case 0:
				val = meanX;
				break;
			case 1:
				val = meanY;
				break;
			case 2:
				val = meanZ;
				break;
			default:
				break;
		}

		val *= 8 * M_PI * M_PI;
		return val;
	}

	return score;
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