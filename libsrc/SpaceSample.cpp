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
#include <hcsrc/Any.h>
#include <string>
#include <iomanip>
#include "ConfSpace.h"
#include "Fibonacci.h"
#include "ConfAxis.h"
#include "Crystal.h"

using namespace HelenCore;

SpaceSample::SpaceSample(ConfSpace *space)
{
	_space = space;
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
}

void SpaceSample::addAverageParameters(RefinementStrategyPtr strategy)
{
	double step = 0.01;
	double tol = 0.00001;
	int n = _space->axisCount();

	for (size_t i = 0; i < n; i++)
	{
		AnyPtr any = AnyPtr(new Any(&_average[i]));
		_anys.push_back(any);
		strategy->addParameter(&*any, Any::get, Any::set, step, tol);
	}
}

void SpaceSample::addDiagonalParameters(RefinementStrategyPtr strategy)
{
	_anys.clear();
	double step = 0.01;
	double tol = 0.00001;
	int n = _space->axisCount();

	for (size_t i = 0; i < n; i++)
	{
		AnyPtr any = AnyPtr(new Any(&_tensor.ptrs[i][i]));
		_anys.push_back(any);

		strategy->addParameter(&*any, Any::get, Any::set, step, tol);
	}
}

void SpaceSample::addTensorParameters(RefinementStrategyPtr strategy)
{
	_anys.clear();
	double step = 0.02;
	double tol = 0.00001;
	int n = _space->axisCount();

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			AnyPtr any = AnyPtr(new Any(&_tensor.ptrs[i][j]));
			_anys.push_back(any);
			
			strategy->addParameter(&*any, Any::get, Any::set, step, tol);
		}
	}
}

void SpaceSample::generatePoints(CrystalPtr crystal, int m)
{
	if (m <= 0)
	{
		m = crystal->getSampleNum();
	}

	_points.clear();
	int dims = _space->axisCount();
	
	Fibonacci fib;
	fib.hyperVolume(dims, m, 1.0);

	std::vector<SpacePoint> hyperpoints = fib.getHyperpoints();
	_points = hyperpoints;
	srand(0);
	
	for (size_t i = 0; i < _points.size(); i++)
	{
		double point = (rand() / (double)RAND_MAX) - 0.5;
		for (size_t j = 0; j < dims; j++)
		{
			_points[i][j] = 0;
			_points[i][j] = point;
//			_points[i][j] = (rand() / (double)RAND_MAX) - 0.5;
//			std::cout << _points[i][j] << " ";
		}
//		std::cout << std::endl;
	}
	
	for (size_t i = 0; i < _points.size(); i++)
	{
//		multMatrix(_tensor, &_points[i][0]);
		
		for (size_t j = 0; j < dims; j++)
		{
			_points[i][j] *= _tensor.ptrs[j][j];
			_points[i][j] += _average[j];
		}
	}
}

double SpaceSample::getWhackDeviation(int resi, int sample)
{
	int n = _space->axisCount();
	
	if (sample >= _points.size())
	{
		return 0;
	}
	
	double sum = 0;
	for (size_t i = 0; i < n; i++)
	{
		ConfAxis *axis = _space->axis(i);
		double add = axis->getWhackDeviationForResidue(resi);
		
		add *= _points[sample][i];
		sum += add;
	}
	
	return sum;
}

double SpaceSample::getTorsionDeviation(int resi, int sample)
{
	int n = _space->axisCount();
	
	if (sample >= _points.size())
	{
		return 0;
	}
	
	double sum = 0;
	for (size_t i = 0; i < n; i++)
	{
		ConfAxis *axis = _space->axis(i);
		double add = axis->getTorsionDeviationForResidue(resi);
		
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


