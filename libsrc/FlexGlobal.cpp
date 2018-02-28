//
//  FlexGlobal.cpp
//  vagabond
//
//  Created by Helen Ginn on 27/12/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "FlexGlobal.h"
#include "AtomGroup.h"
#include "Atom.h"
#include "Bond.h"
#include "Anisotropicator.h"

FlexGlobal::FlexGlobal()
{
	_targetType = FlexTargetMaximiseIsotropy;
	_targetIsoB = 8;
}

double FlexGlobal::matchOriginalBeeScore()
{
	double score = 0;
	double count = 0;

	for (int i = 0; i < _atomGroup->atomCount(); i++)
	{
		AtomPtr atom = _atomGroup->atom(i);

		if (!atom->isFromPDB())
		{
			continue;
		}

		mat3x3 current = atom->getModel()->getRealSpaceTensor();
		mat3x3 original = atom->getTensor();
		mat3x3 diffMat = mat3x3_subtract_mat3x3(current, original);

		double diff = mat3x3_abs_sum_all(diffMat); 

		if (diff != diff) continue;
		if (!isfinite(diff)) continue;

		score += diff;
		count++;

	}

	return score / count;
}

double FlexGlobal::maximiseIsotropyScore()
{
	_atomGroup->propagateChange();
	_atomGroup->refreshPositions();

	int count = 0;

	std::vector<vec3> allPos;

	for (int i = 0; i < _atomGroup->atomCount(); i++)
	{
		AtomPtr atom = _atomGroup->atom(i);

		if (!atom->getModel()->isBond())
		{
			continue;
		}

		BondPtr bond = ToBondPtr(atom->getModel());
		std::vector<BondSample> samples = bond->getFinalPositions();
		vec3 meanPos = bond->getAbsolutePosition();
		allPos.reserve(allPos.size() + samples.size());

		for (int i = 0; i < samples.size(); i++)
		{
			vec3 diff = vec3_subtract_vec3(samples[i].start, meanPos);
			allPos.push_back(diff);
			count++;
		}
	}

	Anisotropicator aniso;
	aniso.setPoints(allPos);
	double score = 0;
	double actualTarget = _targetIsoB / (M_PI * M_PI * 8);

	for (int j = 0; j < 3; j++)
	{
		vec3 anisoAxis = aniso.getAxis(j);
		double sqlength = vec3_sqlength(anisoAxis);
		double diff = fabs(sqlength - actualTarget);
		score += diff;
	}

	return score;
}

void FlexGlobal::maximiseIsotropy()
{
	double sum = 0;
	double count = 0;

	for (int i = 0; i < _atomGroup->atomCount(); i++)
	{
		AtomPtr atom = _atomGroup->atom(i);

		if (!atom->getModel()->isBond())
		{
			continue;
		}

		BondPtr bond = ToBondPtr(atom->getModel());
		double isoTarget = bond->getMeanSquareDeviation();

		sum += isoTarget;
		count++;
	}

	sum /= count;
	_targetIsoB = sum;
	_targetType = FlexTargetMaximiseIsotropy;
}

double FlexGlobal::score(void *object)
{
	FlexGlobal *flexer = static_cast<FlexGlobal *>(object);

	double score = 0;

	switch (flexer->_targetType)
	{
		case FlexTargetMaximiseIsotropy:
		score = flexer->maximiseIsotropyScore();
		break;

		case FlexTargetMatchOrigBFactor:
		score = flexer->matchOriginalBeeScore();
		break;
	}

	return score;
}
