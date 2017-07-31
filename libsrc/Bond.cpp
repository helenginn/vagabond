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

Bond::Bond(AtomPtr major, AtomPtr minor)
{
	_usingTorsion = false;
	_activated = false;
	_major = major;
	_minor = minor;
	_torsionBasis = make_mat3x3();
	_torsionRadians = 0;
	_torsionBlurFromPrev = 0;
	_bendBlur = 0;
	_torsionBlur = 0;
	_bondLength = 0;
	_changedPos = true;
	_changedSamples = true;
	_lastPosition = make_vec3(0, 0, 0);

	vec3 majorPos = getMajor()->getPosition();
	vec3 minorPos = getMinor()->getPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	_bondLength = vec3_length(difference);
	_bondDirection = difference;

//	std::cout << "Bond length, from " << getMajor()->getAtomName()
// << " to " << getMinor()->getAtomName() << " = " << _bondLength << " Å." << std::endl;

	ModelPtr upModel = getMajor()->getModel();

	if (upModel->getClassName() == "Bond")
	{
		std::static_pointer_cast<Bond>(upModel)->addDownstreamAtom(minor);
	}

}

void Bond::addDownstreamAtom(AtomPtr atom)
{
	_downstreamAtoms.push_back(atom);
	vec3 pos = atom->getPosition();
	vec3 start = getMinor()->getPosition();
	vec3 diff = vec3_subtract_vec3(pos, start);

	double angle = vec3_angle_with_vec3(_bondDirection, diff);
	angle -= M_PI / 2;
	double ratio = sin(angle);

	_downRatios.push_back(ratio);
//	std::cout << "Ratio: " << ratio << std::endl;
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

	makeTorsionBasis(hPos, maPos, miPos, lPos, &_torsionRadians);

//	std::cout << "Torsion, on bond " << getMajor()->getAtomName() << " to " <<
//	getMinor()->getAtomName() << " is " << rad2deg(_torsionRadians) << "º" << std::endl;

	_usingTorsion = true;
}

double Bond::getVoxelValue(void *obj, double x, double y, double z)
{
	double distSq = x * x + y * y + z * z;

	if (sqrt(distSq) < 0.16)
	{
		return 1;
	}

	return 0;
}

std::string Bond::getPDBContribution()
{
	std::vector<BondSample> positions = getManyPositions();
	AtomPtr atom = getMinor();
	std::string atomName = atom->getAtomName();
	ElementPtr element = atom->getElement();
	std::string residueName = "LYS";

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 placement = positions[i].start;
		double occupancy = positions[i].occupancy;

		std::cout << "ATOM  ";
		std::cout << "  500 ";
		std::cout << std::setfill(' ') << std::setw(4) << atomName;
		std::cout << " ";
		std::cout << std::setw(3) << residueName;
		std::cout << " A";
		std::cout << "  71";
		std::cout << "    ";
		std::cout << std::fixed << std::setw(8) << std::setprecision(3) << placement.x;
		std::cout << std::setw(8) << std::setprecision(3) << placement.y;
		std::cout << std::setw(8) << std::setprecision(3) << placement.z;
		std::cout << std::setw(6) << std::setprecision(2) << occupancy;
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

	double realScale = 1 / (scale * n);

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 placement = positions[i].start;
		vec3 relative = vec3_subtract_vec3(absolute, placement);
		double occupancy = positions[i].occupancy;

		vec3_mult(&relative, 1 / realScale);

		fft->addToReal(relative.x, relative.y, relative.z, occupancy);
	}

	fft->createFFTWplan(1, false);
	fft->fft(1);
	fft->invertScale();

	FFTPtr fftAbs = _absInherit->getDistribution();
	FFT::multiply(fftAbs, fft);
	FFTPtr check = std::make_shared<FFT>(*fftAbs);
	check->fft(1);

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

	return final;
}

std::vector<BondSample> Bond::sampleMyAngles(double angle, double sigma,
											 double interval)
{
	double sum = 0;
	sigma = fabs(sigma);

	std::vector<BondSample> samples;

	if (sigma <= 0)
	{
		BondSample sample;
		sample.occupancy = 1;
		sample.torsion = angle;
		samples.push_back(sample);

		return samples;
	}

	if (sigma <= interval)
	{
		interval = sigma;
	}

	samples.reserve(sigma * 2 / interval + 1);

	for (double ang = -2 * sigma; ang <= 2 * sigma; ang += interval)
	{
		double relFreq = normal_distribution(ang, sigma);
		sum += relFreq;
		BondSample sample;
		sample.occupancy = relFreq;
		sample.torsion = ang + angle;
		samples.push_back(sample);
	}

	for (int i = 0; i < samples.size(); i++)
	{
		samples[i].occupancy /= sum;
	}

	return samples;
}

std::vector<BondSample> Bond::getCorrelatedAngles(BondSample prev,
												  double lastTorsion,
												  double angle, double blur)
{
	double addBlur = (prev.torsion - lastTorsion) * _torsionBlurFromPrev;
	double newBlur = _torsionRadians + addBlur;

	std::vector<BondSample> set = sampleMyAngles(newBlur, blur);

	return set;
}

std::vector<BondSample> Bond::getManyPositions(bool staticAtom)
{
	if (!_changedSamples && !staticAtom)
	{
//		return _lastSamples;
	}

	vec3 majorPos = getMajor()->getPosition();
	ModelPtr model = getMajor()->getModel();

	std::vector<BondSample> samples;
	
	if (model->getClassName() != "Bond")
	{
		/* We must be connected to something else, oh well */
		/* Torsion basis must be the same. */

		std::vector<BondSample> newSamples;

		double spread = staticAtom ? 0 : _torsionBlur;
		std::vector<BondSample> torsionsOnly = sampleMyAngles(_torsionRadians,
															  spread);
		spread = staticAtom ? 0 : _bendBlur;
		std::vector<BondSample> myBendings = sampleMyAngles(0, _bendBlur,
															ANGLE_SAMPLING);


		vec3 heavyPos = getHeavyAlign()->getPosition();
		vec3 majorPos = getMajor()->getPosition();
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
	std::vector<BondSample> prevSamples = prevBond->getManyPositions(staticAtom);
	AtomPtr prevMajor = prevBond->getMajor();
	double meanLastTorsion = prevBond->getTorsion();

	std::vector<BondSample> newSamples, myTorsions, myBendings;

	double spread = staticAtom ? 0 : _bendBlur;
	myBendings = sampleMyAngles(0, spread, ANGLE_SAMPLING);
	double torsionNumber = prevBond->downstreamAtomNum(getMinor());
	double totalAtoms = prevBond->downstreamAtomCount();

	for (int i = 0; i < prevSamples.size(); i++)
	{
		spread = staticAtom ? 0 : _torsionBlur;
		myTorsions = getCorrelatedAngles(prevSamples[i], meanLastTorsion,
										 _torsionRadians, spread);

		double torsionAngle = prevSamples[i].torsion;
		mat3x3 oldBasis = prevSamples[i].basis;

		double ratio = prevBond->getGeomRatio(torsionNumber);

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

	_lastPosition = getManyPositions(true)[0].start;
	_changedPos = false;
	return _lastPosition;
}

bool Bond::isNotJustForHydrogens()
{
	if (getMinor()->getElement()->electronCount() > 1)
	{
		return true;
	}

	for (int i = 0; i < downstreamAtomCount(); i++)
	{
		AtomPtr atom = downstreamAtom(i);
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

	for (int i = 0; i < downstreamAtomCount(); i++)
	{
		BondPtr bond;
		bond = std::static_pointer_cast<Bond>(downstreamAtom(i)->getModel());
		bond->propagateChange();
	}
}

std::string Bond::description()
{
	std::ostringstream stream;
	stream << "Bond: connects " << getMajor()->getAtomName() << " to "
	<< getMinor()->getAtomName() << std::endl;
	stream << "Bond length: " << _bondLength << " Å" << std::endl;
	stream << "Bond torsion angle: " << rad2deg(_torsionRadians) << std::endl;
	stream << "Bond downstream atoms (" << downstreamAtomCount() << "):" << std::endl;

	for (int i = 0; i < downstreamAtomCount(); i++)
	{
		stream << "\t" << downstreamAtom(i)->getAtomName() << std::endl;
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

	BondPtr newBond = std::static_pointer_cast<Bond>(model);
	int i = newBond->downstreamAtomNum(bond->getMinor());

	if (i >= 0)
	{
		double ratio = newBond->getGeomRatio(i);
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

	BondPtr newBond = std::static_pointer_cast<Bond>(model);
	int i = newBond->downstreamAtomNum(bond->getMinor());

	if (i >= 0)
	{
		double ratio = sin(value - M_PI / 2);
		newBond->setGeomRatio(i, ratio);
	}
}