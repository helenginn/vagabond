//
//  NelderMead.cpp
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "RefinementNelderMead.h"
#include <algorithm>
#include <iostream>

bool NelderMead::converged()
{
	for (int i = 0; i < parameterCount(); i++)
	{
		double limit = _params[i].other_value;
		
		if (_stepMap.count(i) == 0)
		{
			return false;	
		}
		
		double lastDiff = _stepMap[i];
		
		if (fabs(lastDiff) > fabs(limit))
		{
			_lastTag = _params[i].tag;
			return false;	
		}
	}
	
	return true;
}

void NelderMead::addPoints(std::vector<double> *point, std::vector<double> pointToAdd)
{
	for (int i = 0; i < point->size(); i++)
	{
		(*point)[i] += pointToAdd[i];
	}
}

void NelderMead::subtractPoints(std::vector<double> *point, std::vector<double> pointToSubtract)
{
	for (int i = 0; i < point->size(); i++)
	{
		(*point)[i] -= pointToSubtract[i];
	}
}

void NelderMead::scalePoint(std::vector<double> *point, double scale)
{
	for (int i = 0; i < point->size(); i++)
	{
		(*point)[i] *= scale;
	}
}

void NelderMead::setWorstTestPoint(TestPoint &newPoint)
{
	testPoints[testPoints.size() - 1] = newPoint;
}

TestPoint *NelderMead::worstTestPoint()
{
	orderTestPoints();
	return &testPoints[testPoints.size() - 1];
}

std::vector<double> NelderMead::calculateCentroid()
{
	std::vector<double> centroid;
	centroid.resize(parameterCount());
	orderTestPoints();

	for (int i = 0; i < parameterCount(); i++)
	{
		double total = 0;

		for (int j = 0; j < testPoints.size() - 1; j++)
		{
			total += testPoints[j].first[i];
		}

		total /= testPoints.size() - 1;

		centroid[i] = total;
	}

	return centroid;
}

TestPoint NelderMead::reflectOrExpand(std::vector<double> centroid, double scale)
{
	TestPoint *maxPoint = worstTestPoint();

	std::vector<double> diffVec = centroid;
	subtractPoints(&diffVec, maxPoint->first);
	scalePoint(&diffVec, scale);
	std::vector<double> reflectedVec = centroid;
	addPoints(&reflectedVec, diffVec);
	
	for (int i = 0; i < parameterCount(); i++)
	{
		_stepMap[i] = fabs(diffVec[i]);
	}

	TestPoint reflection = std::make_pair(reflectedVec, 0);
	evaluateTestPoint(&reflection);

	return reflection;
}

TestPoint NelderMead::reflectedPoint(std::vector<double> centroid)
{
	return reflectOrExpand(centroid, alpha);
}

TestPoint NelderMead::expandedPoint(std::vector<double> centroid)
{
	return reflectOrExpand(centroid, gamma);
}

TestPoint NelderMead::contractedPoint(std::vector<double> centroid)
{
	return reflectOrExpand(centroid, rho);
}

void NelderMead::reduction()
{
	TestPoint bestPoint = testPoints[0];

	for (int i = 1; i < testPoints.size(); i++)
	{
		TestPoint point = testPoints[i];

		std::vector<double> diffVec = point.first;
		subtractPoints(&diffVec, bestPoint.first);
		scalePoint(&diffVec, sigma);
		std::vector<double> finalVec = bestPoint.first;
		addPoints(&finalVec, diffVec);

		TestPoint contractPoint = std::make_pair(finalVec, 0);
		evaluateTestPoint(&contractPoint);

		testPoints[i] = contractPoint;
	}
}

static bool testPointWorseThanTestPoint(TestPoint one, TestPoint two)
{
	return one.second < two.second;
}

void NelderMead::orderTestPoints()
{
	std::sort(testPoints.begin(), testPoints.end(), testPointWorseThanTestPoint);
}

void NelderMead::evaluateTestPoint(int num)
{
	evaluateTestPoint(&testPoints[num]);
}

void NelderMead::evaluateTestPoint(TestPoint *testPoint)
{
	setTestPointParameters(testPoint);
	double eval = evaluationFunction(evaluateObject);
	testPoint->second = eval;
}

void NelderMead::setTestPointParameters(TestPoint *testPoint)
{
	for (int i = 0; i < parameterCount(); i++)
	{
		setValueForParam(i, testPoint->first[i]);
	}
}

void NelderMead::clearParameters()
{
	RefinementStrategy::clearParameters();

	testPoints.clear();
}

void NelderMead::refine()
{
	RefinementStrategy::refine();

	int testPointCount = (int)parameterCount() + 1;
	testPoints.resize(testPointCount);

	if (parameterCount() == 0)
	return;

	/* Each test point is a vertex? */
	for (int i = 0; i < testPoints.size(); i++)
	{
		testPoints[i].second = 0;
		testPoints[i].first.resize(parameterCount());

		for (int j = 0; j < parameterCount(); j++)
		{
			/* First test point is in the centre */
			if (i == 0)
			{
				testPoints[i].first[j] = getValueForParam(j);
			}

			/* All other test points increase the step size by a
 			 * certain amount in one direction. */
			if (i > 0)
			{
				int minJ = i - 1;
				double scale = 1;

				testPoints[i].first[j] = testPoints[0].first[j] + 
				(j == minJ) * scale * _params[j].step_size;
			}
		}
	}

	for (int i = 0; i < testPoints.size(); i++)
	{
		evaluateTestPoint(i);
	}

	int count = 0;

	while ((!converged() && count < maxCycles))
	{
		std::vector<double> centroid = calculateCentroid();
		count++;

		reportProgress(testPoints[0].second);

		TestPoint reflected = reflectedPoint(centroid);

		if (reflected.second < testPoints[1].second)
		{
			setWorstTestPoint(reflected);
			continue;
		}

		if (reflected.second < testPoints[0].second)
		{
			TestPoint expanded = expandedPoint(centroid);
			bool expandedBetter = (expanded.second < reflected.second);
			setWorstTestPoint(expandedBetter ? expanded : reflected);

			continue;
		}

		TestPoint contracted = contractedPoint(centroid);
		TestPoint *worstPoint = worstTestPoint();

		if (contracted.second < worstPoint->second)
		{
			setWorstTestPoint(contracted);
			continue;
		}
		else
		{
			reduction();
		}
	}

	orderTestPoints();
	reportProgress(testPoints[0].second);
	setTestPointParameters(&testPoints[0]);

	finish();
}

void NelderMead::init()
{
	alpha = 1;
	gamma = 2;
	rho = -0.5;
	sigma = 0.5;
}


