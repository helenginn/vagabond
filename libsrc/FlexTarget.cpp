//
//  FlexTarget.cpp
//  vagabond
//
//  Created by Helen Ginn on 27/12/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "FlexTarget.h"
#include "AtomGroup.h"
#include "Atom.h"
#include "Bond.h"
#include "Anisotropicator.h"

FlexTarget::FlexTarget()
{
	_targetIsoB = 8;
}

double FlexTarget::notStaticScore()
{
	_atomGroup->propagateChange();

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
		std::cout << sqlength << " " << actualTarget << std::endl;
		score += diff;
	}

	std::cout << " = " << score << std::endl;

	return score;
}

double FlexTarget::score(void *object)
{
	return static_cast<FlexTarget *>(object)->notStaticScore();
}
