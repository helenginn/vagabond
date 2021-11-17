// vagabond
// Copyright (C) 2019 Helen Ginn
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

#include "SpaceSample.h"
#include "Anchor.h"
#include "Motion.h"
#include <hcsrc/Any.h>
#include <string>
#include <iomanip>
#include "ConfSpace.h"
#include "Fibonacci.h"
#include "ConfAxis.h"
#include "Superpose.h"
#include "Crystal.h"

using namespace HelenCore;

SpaceSample::SpaceSample(ConfSpace *space)
{
	_used = 3;
	_space = space;
	_superpose = new Superpose();
	int n = _space->axisCount();
	setupMatrix(&_tensor, n, n);
	
	_average = new double[n];
	_stdev = new double[n];
	memset(_average, '\0', sizeof(double) * n);
	memset(_stdev, '\0', sizeof(double) * n);
	
	for (size_t i = 0; i < n; i++)
	{
		_tensor.ptrs[i][i] = 0.;
	}
	
	_average[0] = 0.0;
}

void SpaceSample::addAverageParameters(RefinementStrategyPtr strategy)
{
	double step = 1;
	double tol = 0.0001;
	int n = _space->axisCount();

	for (size_t i = 0; i < n; i++)
	{
		AnyPtr any = AnyPtr(new Any(&_average[i], 0.005));
		_anys.push_back(any);
		strategy->addParameter(&*any, Any::get, Any::set, step, tol);
	}
}

void SpaceSample::addDiagonalParameters(RefinementStrategyPtr strategy)
{
	_anys.clear();
	double step = 1;
	double tol = 0.0001;
	int n = _space->axisCount();

	for (size_t i = 0; i < n; i++)
	{
		AnyPtr any = AnyPtr(new Any(&_tensor.ptrs[i][i], 0.005));
		_anys.push_back(any);

		strategy->addParameter(&*any, Any::get, Any::set, step, tol);
	}
}

void SpaceSample::fillInTensorGaps()
{
	int dims = _space->axisCount();

	for (size_t i = 0; i < dims; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			_tensor.ptrs[i][j] = _tensor.ptrs[j][i];
		}
	}
}

void SpaceSample::savePositions()
{
	superpose()->savePositions();
}

void SpaceSample::addTensorParameters(RefinementStrategyPtr strategy)
{
	_anys.clear();
	double step = 1;
	double tol = 0.0001;
	int n = _space->axisCount();

	for (size_t i = 0; i < n && i < _used; i++)
	{
		for (size_t j = i; j < n; j++)
		{
			AnyPtr any = AnyPtr(new Any(&_tensor.ptrs[i][j], 0.005));
			_anys.push_back(any);
			
			strategy->addParameter(&*any, Any::get, Any::set, step, tol);
		}
	}
}

void SpaceSample::generatePoints(TotalModelPtr total, bool use)
{
	int m = total->getSampleNum() * 1.2;

	_points.clear();
	int dims = _space->axisCount();
	
	Fibonacci fib;
	fib.hyperVolume(dims, _used, m, 1.0);

	std::vector<SpacePoint> hyperpoints = fib.getHyperpoints();
	_points = hyperpoints;
	_original = hyperpoints;
	
	for (size_t i = 0; i < motionCount(); i++)
	{
		double scale = _anchors[i]->getRMSD() * 2;
		getMotion(i)->prepareMotions(hyperpoints, scale);
	}

	srand(0);
	
	if (!use)
	{
		return;
	}
	
	fillInTensorGaps();
	
	for (size_t i = 0; i < _points.size(); i++)
	{
		multMatrix(_tensor, &_points[i][0]);
		
		for (size_t j = 0; j < dims; j++)
		{
			_points[i][j] += _average[j];
		}
	}
	
	for (size_t i = 0; i < motionCount(); i++)
	{
	}
}

double SpaceSample::getDeviation(int resi, int sample, bool phi)
{
	int n = _space->axisCount();
	
	if (sample >= _points.size())
	{
		return 0;
	}
	
	return _space->weightedTorsion(resi, phi, _points[sample]);
	
	double sum = 0;
	for (size_t i = 0; i < n; i++)
	{
		ConfAxis *axis = _space->axis(i);
		double add = 0;
		if (phi)
		{
			add = axis->getPhiDeviationForResidue(resi);
		}
		else
		{
			add = axis->getPsiDeviationForResidue(resi);
		}
		
		add *= _points[sample][i];
		sum += add;
	}
	
	return sum;
}

SpaceSample::~SpaceSample()
{
	delete [] _average;
	delete [] _stdev;

	freeMatrix(&_tensor);
}

vec3 SpaceSample::point3D(int i)
{
	vec3 point = empty_vec3();
	
	for (size_t j = 0; j < 3; j++)
	{
		*(&point.x + j) = _original[i][j];
	}
	
	return point;
}

void SpaceSample::calculateSuperpositions()
{
	superpose()->calculateDeviations();
}

void SpaceSample::setAtoms(AtomList atoms)
{
	superpose()->setAtoms(atoms);
}
