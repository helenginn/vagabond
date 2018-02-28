//
//  RefinementSnake.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/09/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "RefinementSnake.h"
#include "Bond.h"
#include "Sampler.h"
#include <float.h>
#include <sstream>

std::vector<double> RefinementSnake::getTorsionList()
{
	std::vector<double> oldTorsions;

	for (int i = 0; i < bondCount(); i++)
	{
		oldTorsions.push_back(Bond::getTorsion(&*bond(i)));
	}

	return oldTorsions;
}


void RefinementSnake::applyTorsionList(std::vector<double> list)
{
	for (int i = 0; i < bondCount(); i++)
	{
		Bond::setTorsion(&*bond(i), list[i]);
	}
}

void RefinementSnake::setParentSampler(Sampler *sampler)
{
	_parentSampler = sampler;
}

void RefinementSnake::refine()
{
	Sampler sampler;
	std::vector<double> oldTorsions = getTorsionList();

	double startValue = Sampler::score(_parentSampler);
	double nextValue = FLT_MAX;

	std::vector<int> nTries(bondCount(), 0);
	int tryPointer = bondCount() - 1;
	std::ostringstream results;

	for (int i = 0; i < bondCount(); i++)
	{
		bond(i)->copyTarget(_parentSampler);
		bond(i)->setupSampling();
		double value = bond(i)->sample(false);
		results << value << ", ";
	}

	nextValue = Sampler::score(_parentSampler);

	while (nextValue > startValue)
	{
		results.str("");
		results.clear();

		if (nTries[tryPointer] < 3)
		{
			nTries[tryPointer]++;
			std::vector<double> nextList;
			nextList = bond(tryPointer)->getNextResult(nTries[tryPointer]);

			//		std::cout << "Trying next attempt... (" << tryPointer
			//		<< "/" << nTries[tryPointer] << ")" << std::endl;

			Bond::setTorsion(&*bond(tryPointer), nextList[0]);

			for (int i = tryPointer + 1; i < nTries.size(); i++)
			{
				nTries[i] = 0;
				bond(i)->copyTarget(_parentSampler);
				bond(i)->setupSampling();
				double value = bond(i)->sample(false);
			}
		}
		else
		{
			tryPointer--;

			if (tryPointer < 0)
			{
				break;
			}
		}

		nextValue = Sampler::score(_parentSampler);
	}

	if (nextValue > startValue)
	{
		applyTorsionList(oldTorsions);

		if (!_silent)
		{
			std::cout << "No change for bond " << _parentSampler->getJobName() << " (" << startValue << ") ";
			std::cout << /*"debug " << results.str() <<*/ std::endl;
		}
	}
	else
	{
		if (!_silent)
		{
			std::cout << "Snake change for bond " << _parentSampler->getJobName() << " (" << startValue << " to " << nextValue
			<< ")" << std::endl;
		}
	}
}
