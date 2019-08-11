// Vagabond : bond-based macromolecular model refinement
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

#ifndef __vagabond__NelderMead__
#define __vagabond__NelderMead__

#include <stdio.h>
#include "RefinementStrategy.h"
#include <map>

typedef std::pair<std::vector<double>, double> TestPoint;
typedef std::map<int, double> StepMap;

class RefinementNelderMead : public RefinementStrategy
{
public:
	void init();
	RefinementNelderMead() : RefinementStrategy() { init(); };
	virtual void refine();

	virtual void clearParameters();
private:
	double alpha;
	double gamma;
	double rho;
	double sigma;
	
	StepMap _stepMap;
	std::string _lastTag;

	std::vector<TestPoint> testPoints;

	void setWorstTestPoint(TestPoint &newPoint);
	TestPoint *worstTestPoint();
	void orderTestPoints();
	void evaluateTestPoint(int num);
	void evaluateTestPoint(TestPoint *testPoint);
	void setTestPointParameters(TestPoint *testPoint);
	std::vector<double> calculateCentroid();

	TestPoint reflectOrExpand(std::vector<double> centroid, double scale);
	TestPoint reflectedPoint(std::vector<double> centroid);
	TestPoint expandedPoint(std::vector<double> centroid);
	TestPoint contractedPoint(std::vector<double> centroid);
	void reduction();
	bool converged();

	void addPoints(std::vector<double> *point, std::vector<double> pointToAdd);
	void scalePoint(std::vector<double> *point, double scale);
	void subtractPoints(std::vector<double> *point, std::vector<double> pointToSubtract);
};

#endif /* defined(__vagabond__NelderMead__) */
