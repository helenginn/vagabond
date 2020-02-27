// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "RefinementList.h"
#include <iostream>
#include <iomanip>

RefinementList::RefinementList() : RefinementStrategy()
{
	_cycleNum = 0;

}

void RefinementList::addTestSet(std::vector<double> &vals)
{
	if (parameterCount() != vals.size())
	{
		std::cout << "Test set has different parameter number" << std::endl;
		return;
	}

	_tests.push_back(vals);
}

void RefinementList::applyTest(int num)
{
	std::vector<double> test = _tests[num];

	for (int i = 0; i < test.size(); i++)
	{
		setValueForParam(i, test[i]);
	}
}

void RefinementList::refine()
{
	RefinementStrategy::refine();
	int bestCycle = 0;
	double bestScore = _prevScore;
	
	while (_cycleNum < _tests.size())
	{
		applyTest(_cycleNum);

		double eval = (*evaluationFunction)(evaluateObject);

		if (eval < bestScore)
		{
			bestScore = eval;
			bestCycle = _cycleNum;
		}
		reportProgress(bestScore);

		_cycleNum++;
	}

	_chosen = bestCycle;

	applyTest(bestCycle);
	finish();
}
