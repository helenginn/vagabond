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

Bond::Bond(AtomPtr major, AtomPtr minor)
{
	_usingTorsion = false;
	_activated = false;
	_major = major;
	_minor = minor;
	_torsionBasis = make_mat3x3();
	_torsionRadians = 0;
	_bondLength = 0;
	_minorGeometry = BondGeometryTetrahedral;

	vec3 majorPos = getMajor()->getPosition();
	vec3 minorPos = getMinor()->getPosition();

	vec3 difference = vec3_subtract_vec3(majorPos, minorPos);
	_bondLength = vec3_length(difference);
	_bondDirection = difference;

	std::cout << "Bond length: " << _bondLength << " Å." << std::endl;

	ModelPtr upModel = getMajor()->getModel();

	if (upModel->getClassName() == "Bond")
	{
		upModel->addDownstreamAtom(minor);
	}

}

void Bond::activate(AtomGroupPtr group)
{
	getMinor()->setModel(shared_from_this());

	if (group)
	{
		BondPtr myself = std::static_pointer_cast<Bond>(shared_from_this());
		group->addBond(myself);
	}
	_activated = true;
}

mat3x3 Bond::makeTorsionBasis(vec3 _specificDirection, double *newAngle)
{
	vec3 pPos = _heavyAlign.lock()->getPosition();
	vec3 aPos = getMajor()->getPosition();
	vec3 bPos = getMinor()->getPosition();
	vec3 lPos = _lightAlign.lock()->getPosition();

	vec3 a2p = vec3_subtract_vec3(pPos, aPos);
	vec3 a2l = vec3_subtract_vec3(lPos, aPos);
	vec3 a2b = vec3_subtract_vec3(bPos, aPos);

	double sql = vec3_sqlength(a2b);
	double dotp = vec3_dot_vec3(a2p, a2b);
	double dotl = vec3_dot_vec3(a2l, a2b);
	double distp = dotp / sql;
	double distl = dotl / sql;

	vec3 a2bp = a2b;
	vec3 a2bl = a2b;
	vec3_mult(&a2bp, distp);
	vec3_mult(&a2bl, distl);
	vec3 heavy_join = vec3_add_vec3(aPos, a2bp);
	vec3 light_join = vec3_add_vec3(aPos, a2bl);

	vec3 reverse_bond = _bondDirection;
	vec3_mult(&reverse_bond, -1);
	vec3 xNew = vec3_subtract_vec3(pPos, heavy_join);
	vec3 light = vec3_subtract_vec3(lPos, light_join);
	mat3x3 test = mat3x3_rhbasis(xNew, reverse_bond);
	_torsionBasis = test;

	double angle = vec3_angle_with_vec3(light, xNew);

	vec3 recross = vec3_cross_vec3(light, xNew);
	double cosine = vec3_cosine_with_vec3(recross, reverse_bond);

	if (cosine > 0)
	{
		angle += M_PI;
	}

	std::cout << "Torsion angle: " << angle * 180 / M_PI << "º."<< std::endl;
	_torsionRadians = angle;

	return test;
}

void Bond::setTorsionAtoms(AtomPtr heavyAlign, AtomPtr lightAlign)
{
	_heavyAlign = heavyAlign;
	_lightAlign = lightAlign;

	/* Make torsion basis.
	 * Make any starting set of angles with correct Z axis. */

	makeTorsionBasis(_bondDirection);

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
	/* Getting the distribution for the 'major' atom */

	FFTPtr inherited = getMajor()->getBlur();

	double n = ATOM_SAMPLING_COUNT;
	double scale = 1 / (2.0 * MAX_SCATTERING_DSTAR);

	/* Deal with a needed width of 4.0 Å for now */

	prepareDistribution(n, scale, this, &getVoxelValue);
	FFTPtr mine = getDistributionCopy();
	mine->fft(1);
	FFT::multiply(mine, inherited);

	return mine;
}

/* This gets the position of the minor atom. */
vec3 Bond::getPosition()
{
	/* You'll need to have some "if changed..." then reconstruct things. */

	vec3 majorPos = getMajor()->getPosition(); // my starting atom position!
	ModelPtr model = getMajor()->getModel();
	double myLength = getBondLength(this);     // my bond length!

	if (model->getClassName() != "Bond")
	{
		/* We must be connected to something else, oh well */

		/* Fix me and include sensible calculation! */

		return vec3_subtract_vec3(majorPos, _bondDirection);
	}

	BondPtr prevBond = std::static_pointer_cast<Bond>(model);

	double torsionAngle = prevBond->getTorsion(); // my torsion angle!
	mat3x3 torsionBasis = prevBond->getTorsionBasis(); // my torsion vectors!
	double ratio = 1./3.; // replace me. My x/z ratio!
	double torsionNumber = prevBond->downstreamAtomNum(getMinor());
	double totalAtoms = prevBond->downstreamAtomCount();

	if (torsionNumber < 0)
	{
		shout_at_helen("Something has gone horrendously wrong\n"\
					   "in the calculation of torsion angle.");
	}

	torsionAngle += deg2rad(360) * torsionNumber / totalAtoms;

	/* Calculate the right ratio of x-to-z from the major atom. */
	vec3 atomWrtMajor = make_vec3(1, 0, ratio);

	/* Set these to adhere to the correct bond length. */
	vec3_set_length(&atomWrtMajor, myLength);

	/* Rotate the bond so it lines up with the torsion angle. */
	vec3 zAxis = make_vec3(0, 0, 1);
	mat3x3 torsion_turn = mat3x3_unit_vec_rotation(zAxis, -torsionAngle);
	mat3x3_mult_vec(torsion_turn, &atomWrtMajor);

	/* Reset the basis vectors so that they line up with the previous
	 * bond and the torsion angle is correct. */
	mat3x3_mult_vec(torsionBasis, &atomWrtMajor);

	/* Add this to the major atom position! */
	vec3 final = vec3_add_vec3(majorPos, atomWrtMajor);

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