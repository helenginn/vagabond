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

Bond::Bond(AtomPtr major, AtomPtr minor)
{
	_usingTorsion = false;
	_activated = false;
	_major = major;
	_minor = minor;
	_torsionBasis = make_mat3x3();
	_torsionRadians = 0;
	_torsionBlur = 0;
	_bondLength = 0;
	_minorGeometry = BondGeometryTetrahedral;

	vec3 majorPos = getMajor()->getPosition();
	vec3 minorPos = getMinor()->getPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	_bondLength = vec3_length(difference);
	_bondDirection = difference;

	std::cout << "Bond length = " << _bondLength << " Å." << std::endl;

	ModelPtr upModel = getMajor()->getModel();

	if (upModel->getClassName() == "Bond")
	{
		upModel->addDownstreamAtom(minor);
	}

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

mat3x3 Bond::makeTorsionBasis(vec3 _specificDirection, vec3 hPos, vec3 maPos,
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

	vec3 reverse_bond = _bondDirection;
	vec3_mult(&reverse_bond, -1);
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
			angle += M_PI;
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

	makeTorsionBasis(_bondDirection, hPos, maPos,
					 miPos, lPos, &_torsionRadians);

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

FFTPtr Bond::getDistribution()
{
	std::vector<BondSample> positions = getManyPositions();
	vec3 absolute = getMinor()->getPosition();

	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2.0 * MAX_SCATTERING_DSTAR);

	FFTPtr fft = FFTPtr(new FFT());
	fft->create(n);
	fft->setScales(scale);

	for (int i = 0; i < positions.size(); i++)
	{
		vec3 placement = positions[i].start;
		vec3 relative = vec3_subtract_vec3(absolute, placement);
		double occupancy = positions[i].occupancy;

		vec3_mult(&relative, 1 / scale);

		fft->addToReal(relative.x, relative.y, relative.z, occupancy);
	}

	fft->createFFTWplan(1, false);
	fft->fft(1);

	FFTPtr fftAbs = _absInherit->getDistribution();
	FFT::multiply(fftAbs, fft);

	return fftAbs;

	/* Getting the distribution for the 'major' atom */

	FFTPtr inherited = getMajor()->getBlur();


	/* Deal with a needed width of 4.0 Å for now */

	prepareDistribution(n, scale, this, &getVoxelValue);
	FFTPtr mine = getDistributionCopy();
	mine->fft(1);
	FFT::multiply(mine, inherited);

	return mine;
}

vec3 Bond::positionFromTorsion(mat3x3 torsionBasis, double angle,
							   double ratio, vec3 start, vec3 *heavyPtr,
							   mat3x3 *myBasis)
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

	if (myBasis && heavyPtr)
	{
		vec3 heavy = *heavyPtr;
		vec3 none = make_vec3(0, 0, 0);

		*myBasis = makeTorsionBasis(atomWrtMajor, heavy, start, final,
										   none);
	}

	return final;
}

std::vector<BondSample> Bond::sampleMyPositions()
{
	double interval = deg2rad(0.5);
	double sum = 0;
	double sigma = fabs(_torsionBlur);

	std::vector<BondSample> samples;

	if (_torsionBlur <= interval)
	{
		BondSample sample;
		sample.occupancy = 1;
		sample.torsion = _torsionRadians;
		samples.push_back(sample);

		return samples;
	}

	for (double ang = -sigma; ang <= sigma; ang += interval)
	{
		double relFreq = normal_distribution(ang, 0, sigma);
		sum += relFreq;
		BondSample sample;
		sample.occupancy = relFreq;
		sample.torsion = ang + _torsionRadians;
		samples.push_back(sample);
	}

	for (int i = 0; i < samples.size(); i++)
	{
		samples[i].occupancy /= sum;
	}

	return samples;
}

std::vector<BondSample> Bond::getManyPositions()
{
	vec3 majorPos = getMajor()->getPosition();
	ModelPtr model = getMajor()->getModel();

	std::vector<BondSample> samples;
	
	if (model->getClassName() != "Bond")
	{
		/* We must be connected to something else, oh well */
		/* Torsion basis must be the same. */

		std::vector<BondSample> torsionsOnly = sampleMyPositions();
		vec3 majorPos = getMajor()->getPosition();
		vec3 start = vec3_subtract_vec3(majorPos, _bondDirection);

		for (int i = 0; i < torsionsOnly.size(); i++)
		{
			torsionsOnly[i].basis = _torsionBasis;
			torsionsOnly[i].start = start;
		}

		return torsionsOnly;
	}

	BondPtr prevBond = std::static_pointer_cast<Bond>(model);
	std::vector<BondSample> prevSamples = prevBond->getManyPositions();
	double ratio = prevBond->getAtomicAngle();
	AtomPtr prevMajor = prevBond->getMajor();

	std::vector<BondSample> newSamples, mySamples;
	BondSample aSample;
	aSample.torsion = _torsionRadians;
	aSample.occupancy = 1;
	mySamples.push_back(aSample);

	for (int i = 0; i < prevSamples.size(); i++)
	{
		double torsionAngle = prevSamples[i].torsion;
		mat3x3 oldBasis = prevSamples[i].basis;

		double torsionNumber = prevBond->downstreamAtomNum(getMinor());
		double totalAtoms = prevBond->downstreamAtomCount();

		if (torsionNumber < 0)
		{
			shout_at_helen("Something has gone horrendously wrong\n"\
						   "in the calculation of torsion angle.");
		}

		vec3 heavy = prevMajor->getPosition();
		torsionAngle += deg2rad(360) * torsionNumber / totalAtoms;

		mat3x3 newBasis;
		vec3 final = positionFromTorsion(oldBasis, torsionAngle, ratio,
										 majorPos, &heavy, &newBasis);

		for (int j = 0; j < mySamples.size(); j++)
		{
			BondSample nextSample;
			nextSample.basis = newBasis;
			nextSample.start = final;
			nextSample.torsion = mySamples[j].torsion;
			nextSample.occupancy = mySamples[j].occupancy *
			prevSamples[i].occupancy;

			newSamples.push_back(nextSample);
		}
	}

	return newSamples;
}

/* This gets the position of the minor atom. */
vec3 Bond::getStaticPosition()
{
	vec3 majorPos = getMajor()->getPosition();
	ModelPtr model = getMajor()->getModel();

	if (model->getClassName() != "Bond")
	{
		/* We must be connected to something else, oh well */

		return vec3_subtract_vec3(majorPos, _bondDirection);
	}

	BondPtr prevBond = std::static_pointer_cast<Bond>(model);

	double torsionAngle = prevBond->getTorsion(); // my torsion angle!
	mat3x3 torsionBasis = prevBond->getTorsionBasis(); // my torsion vectors!
	double ratio = prevBond->getAtomicAngle(); // replace me. My x/z ratio!
	double torsionNumber = prevBond->downstreamAtomNum(getMinor());
	double totalAtoms = prevBond->downstreamAtomCount();

	if (torsionNumber < 0)
	{
		shout_at_helen("Something has gone horrendously wrong\n"\
					   "in the calculation of torsion angle.");
	}

	torsionAngle += deg2rad(360) * torsionNumber / totalAtoms;

	vec3 final = positionFromTorsion(torsionBasis, torsionAngle, ratio,
									majorPos);

	return final;
}

bool Bond::isNotJustForHydrogens()
{
	if (getMinor()->getElement()->electronCount() <= 1)
	{
		return false;
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