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

Bond::Bond(AtomPtr major, AtomPtr minor, int group)
{
	_usingTorsion = false;
	_activated = false;
	_major = major;
	_minor = minor;
	_activeGroup = 0;
	_torsionBasis = make_mat3x3();
	_torsionAngles.push_back(0);
	_torsionBlurs.push_back(0);
	_torsionBlurFromPrev = 0;
	_bendBlur = 0;
	_bondLength = 0;
	_changedPos = true;
	_changedSamples = true;
	_lastPosition = make_vec3(0, 0, 0);

	vec3 majorPos = getMajor()->getPosition();
	vec3 minorPos = getMinor()->getPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	_bondLength = vec3_length(difference);

	deriveBondLength();
	vec3_set_length(&difference, _bondLength);

	_bondDirection = difference;

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
	_torsionBasis = other._torsionBasis;
	_torsionAngles = other._torsionAngles;
	_torsionBlurFromPrev = 0;
	_bendBlur = 0;
	_torsionBlurs = std::vector<double>(1, 0);
	_bondLength = other._bondLength;
	_changedPos = true;
	_changedSamples = true;
	_lastPosition = make_vec3(0, 0, 0);
	_bondDirection = other._bondDirection;
	_absInherit = other._absInherit;
	_heavyAlign = other._heavyAlign;
	_lightAlign = other._lightAlign;
	_bendToAtom = other._bendToAtom;
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
		std::cout << "Setting bond length to " << _bondLength << " Å." << std::endl;
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

	AtomType type3 = downstreamAtom(group, n)->getGeomType();

	GeomTable table = GeomTable::getGeomTable();
	double angle = table.getBondAngle(type1, type2, type3);

	if (angle > 0)
	{
		double ratio = sin(angle - M_PI / 2);
		setGeomRatio(group, n, ratio);

		std::cout << "Setting angle to " << rad2deg(angle) << std::endl;
	}
}

void Bond::addDownstreamAtom(AtomPtr atom, int group)
{
	while (_downstreamAtoms.size() <= group)
	{
		_downstreamAtoms.push_back(AtomList());
		_torsionAngles.push_back(0);
		_torsionBlurs.push_back(0);
		_downRatios.push_back(std::vector<double>());
	}

	_downstreamAtoms[group].push_back(atom);

	vec3 pos = atom->getInitialPosition();
	vec3 start = getMinor()->getInitialPosition();
	vec3 diff = vec3_subtract_vec3(pos, start);

	double angle = vec3_angle_with_vec3(_bondDirection, diff);
	angle -= M_PI / 2;
	double ratio = sin(angle);

	_downRatios[group].push_back(ratio);

	deriveBondAngle(_downRatios.size() - 1, _downRatios[group].size() - 1);
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

void Bond::activate(AtomGroupPtr group, AbsolutePtr inherit)
{
	getMinor()->setModel(shared_from_this());

	if (group)
	{
		BondPtr myself = std::static_pointer_cast<Bond>(shared_from_this());
		group->addBond(myself);
	}

	if (inherit)
	{
		setAbsoluteInheritance(inherit);
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
	_torsionBasis = test;

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

	return test;
}

void Bond::setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign)
{
	_heavyAlign = heavyAlign;
	_lightAlign = lightAlign;

	/* Make torsion basis.
	 * Make any starting set of angles with correct Z axis. */

	vec3 hPos = _heavyAlign.lock()->getPosition();
	vec3 maPos = getMajor()->getPosition();
	vec3 miPos = getMinor()->getPosition();
	vec3 lPos = _lightAlign.lock()->getPosition();

	makeTorsionBasis(hPos, maPos, miPos, lPos, &_torsionAngles[0]);

//	std::cout << "Torsion, on bond " << getMajor()->getAtomName() << " to " <<
//	getMinor()->getAtomName() << " is " << rad2deg(_torsionAngle) << "º" << std::endl;

	_usingTorsion = true;
}

std::string Bond::getPDBContribution()
{
	AtomPtr atom = getMinor();
	std::string atomName = atom->getAtomName();
	ElementPtr element = atom->getElement();
	std::string residueName = "LYS";

	if (element->getSymbol() == "H")
	{
		return std::string();
	}

	const int tries = 10;

	for (int i = 0; i < tries; i++)
	{
		std::vector<BondSample> positions = getManyPositions(false, true);

	//	std::cout << "Confirming " << positions.size() << " size" << std::endl;
		vec3 placement = positions[0].start;
		double occupancy = positions[0].occupancy;

		std::cout << "ATOM  ";
		std::cout << "  500";
		std::cout << std::setfill(' ') << std::setw(4) << atomName;
		std::cout << "  ";
		std::cout << std::setw(3) << residueName;
		std::cout << " A";
		std::cout << " 123";
		std::cout << "    ";
		std::cout << std::fixed << std::setw(8) << std::setprecision(3) << placement.x;
		std::cout << std::setw(8) << std::setprecision(3) << placement.y;
		std::cout << std::setw(8) << std::setprecision(3) << placement.z;
		std::cout << std::setw(6) << std::setprecision(2) << occupancy / double(tries);
		std::cout << std::setw(6) << std::setprecision(2) << _absInherit->getBFactor();
		std::cout << "          ";
		std::cout << std::setw(2) << element->getSymbol();
		std::cout << "  " << std::endl;;
	}

	return std::string();
}

FFTPtr Bond::getDistribution()
{
	std::vector<BondSample> positions = getManyPositions();
	vec3 absolute = getMinor()->getPosition();

	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2 * MAX_SCATTERING_DSTAR);

	double realLimits = (scale * n);

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 placement = positions[i].start;
		vec3 relative = vec3_subtract_vec3(placement, absolute);
		double occupancy = positions[i].occupancy;

		vec3_mult(&relative, 1 / realLimits);

		if (getMinor()->getAtomName() == "CG")
		{
		//	std::cout << vec3_desc(relative) << std::endl;
		}

		fft->addToReal(relative.x, relative.y, relative.z, occupancy);
	}

	fft->createFFTWplan(1, false);
	fft->fft(1);
	fft->invertScale();

	FFTPtr fftAbs = _absInherit->getDistribution();
	FFT::multiply(fftAbs, fft);

	return fftAbs;
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

//	std::cout << "Angle: " << angle << " to " << vec3_desc(final) << std::endl;

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

	//	std::cout << "Angle " << rad2deg(angle) << ", sigma " << rad2deg(sigma)
	//	<< " generates " << rad2deg(sample.torsion) << "º." << std::endl;
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
		interval = sigma;
	}

	samples.reserve(sigma * 2 / interval + 1);
	int count = 0;

	for (double ang = -2 * sigma; ang <= 2 * sigma; ang += interval)
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

std::vector<BondSample> Bond::getCorrelatedAngles(BondSample prev,
												  double lastTorsion,
												  double angle, double blur,
												  bool singleState, int group)
{
	double addBlur = (prev.torsion - lastTorsion) * _torsionBlurFromPrev;
	double newBlur = _torsionAngles[group] + addBlur;

	std::vector<BondSample> set = sampleMyAngles(newBlur, blur, singleState);

	return set;
}

std::vector<BondSample> Bond::getManyPositions(bool staticAtom,
											   bool singleState,
											   int group)
{
	if (!_changedSamples && !staticAtom && !singleState)
	{
//		return _lastSamples;
	}

	ModelPtr model = getMajor()->getModel();

	std::vector<BondSample> samples;
	vec3 majorPos = getMajor()->getPosition();

	if (model->getClassName() != "Bond")
	{
		/* We must be connected to something else, oh well */
		/* Torsion basis must be the same. */

		std::vector<BondSample> newSamples;

		double spread = staticAtom ? 0 : _torsionBlurs[group];
		std::vector<BondSample> torsionsOnly, myBendings;
		torsionsOnly = sampleMyAngles(_torsionAngles[group], spread, singleState);
		spread = staticAtom ? 0 : _bendBlur;

		myBendings = sampleMyAngles(0, spread, singleState);


		vec3 heavyPos = getHeavyAlign()->getPosition();
		vec3 none = {0, 0, 1};
		if (!_bendToAtom.lock())
		{
			shout_at_helen("Need to specify which atom you bend to\n"\
						   "if you want to bend the bond.");
		}

		vec3 bendToPos = _bendToAtom.lock()->getPosition();
		vec3 toBendDir = vec3_subtract_vec3(majorPos, bendToPos);

		for (int j = 0; j < myBendings.size(); j++)
		{
			vec3 cross = vec3_cross_vec3(_bondDirection, toBendDir);
			vec3_set_length(&cross, 1);
			mat3x3 bend = mat3x3_unit_vec_rotation(cross, myBendings[j].torsion);
			vec3 newDirection = _bondDirection;
			mat3x3_mult_vec(bend, &newDirection);
			vec3 start = vec3_subtract_vec3(majorPos, newDirection);

			mat3x3 newBasis = makeTorsionBasis(heavyPos, majorPos, start, none);

			for (int i = 0; i < torsionsOnly.size(); i++)
			{
				BondSample newSample;
				newSample.basis = newBasis;
				newSample.start = start;
				newSample.torsion = torsionsOnly[i].torsion;
				newSample.occupancy = torsionsOnly[i].occupancy *
				myBendings[j].occupancy;
				newSamples.push_back(newSample);
			}
		}

	//	_lastSamples = newSamples;
		_changedSamples = false;

		return newSamples;
	}

	BondPtr prevBond = std::static_pointer_cast<Bond>(model);
	int myGroup = -1;
	double torsionNumber = prevBond->downstreamAtomNum(getMinor(), &myGroup);
	double totalAtoms = prevBond->downstreamAtomCount(myGroup);

	std::vector<BondSample> prevSamples = prevBond->getManyPositions(staticAtom,
																	 singleState,
																	 myGroup);

	AtomPtr prevMajor = prevBond->getMajor();
	double meanLastTorsion = getTorsion(&*prevBond);

	std::vector<BondSample> newSamples, myTorsions, myBendings;

	double spread = staticAtom ? 0 : _bendBlur;
	myBendings = sampleMyAngles(0, spread, singleState);

	for (int i = 0; i < prevSamples.size(); i++)
	{
		spread = staticAtom ? 0 : _torsionBlurs[group];
		myTorsions = getCorrelatedAngles(prevSamples[i], meanLastTorsion,
										 _torsionAngles[group], spread,
										 singleState);

		double torsionAngle = prevSamples[i].torsion;
		mat3x3 oldBasis = prevSamples[i].basis;

		double ratio = prevBond->getGeomRatio(myGroup, torsionNumber);

		if (torsionNumber < 0)
		{
			shout_at_helen("Something has gone horrendously wrong\n"\
						   "in the calculation of torsion angle.\n" +
						   description());
		}

		vec3 heavy = prevMajor->getPosition();
		torsionAngle += deg2rad(360) * torsionNumber / totalAtoms;

		mat3x3 newBasis;

		/* Prepping for bending */
		vec3 zAxis = {0, 0, 1};
		mat3x3_mult_vec(oldBasis, &zAxis);
		const vec3 none = {0, 0, 1};

		/* do the bendy thing */

		for (int j = 0; j < myBendings.size(); j++)
		{
			vec3 torsionPos = positionFromTorsion(oldBasis, torsionAngle,
												  ratio, majorPos);

			vec3 cross = vec3_cross_vec3(torsionPos, zAxis);
			vec3_set_length(&cross, 1);
			mat3x3 bend = mat3x3_unit_vec_rotation(cross, myBendings[j].torsion);
			mat3x3_mult_vec(bend, &torsionPos);

			newBasis = makeTorsionBasis(heavy, majorPos, torsionPos, none);

			for (int k = 0; k < myTorsions.size(); k++)
			{
				BondSample nextSample;
				nextSample.basis = newBasis;
				nextSample.start = torsionPos;
				nextSample.torsion = myTorsions[k].torsion;
				nextSample.occupancy = myTorsions[k].occupancy *
				myBendings[j].occupancy * prevSamples[i].occupancy;

				newSamples.push_back(nextSample);
			}
		}
	}

//	_lastSamples = newSamples;
	_changedSamples = false;

	return newSamples;
}

/* This gets the position of the minor atom. */
vec3 Bond::getStaticPosition()
{
	if (!_changedPos)
	{
		return _lastPosition;
	}

	_lastPosition = getManyPositions(true, false)[0].start;
	_changedPos = false;
	return _lastPosition;
}

bool Bond::isNotJustForHydrogens()
{
	if (getMinor()->getElement()->electronCount() > 1)
	{
		return true;
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
	stream << "Bond torsion angle: " << rad2deg(_torsionAngles[0]) << std::endl;
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
		double angle = asin(ratio) + M_PI / 2;
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
		double ratio = sin(value - M_PI / 2);
		newBond->setGeomRatio(myGroup, i, ratio);
	}

	static_cast<Bond *>(object)->propagateChange();
}

bool Bond::splitBond()
{
	BondPtr me = std::static_pointer_cast<Bond>(shared_from_this());
	BondPtr parent = std::static_pointer_cast<Bond>(getParentModel());
	int last = downstreamAtomGroupCount();
	me->duplicateDownstream(parent, last);
	parent->setActiveGroup(last);
	double torsion = getTorsion(&*parent);
	torsion += deg2rad(180);

	if (torsion > deg2rad(360))
	{
		torsion -= deg2rad(360);
	}

	parent->setTorsion(this, torsion);
	parent->setActiveGroup(0);

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

		if (model->getClassName() != "Bond")
		{
			continue;
		}

		nextBond->duplicateDownstream(duplBond, 0);
	}
}