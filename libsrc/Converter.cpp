// Vagabond
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

#include "Converter.h"
#include "../libica/svdcmp.h"
#include <cstring>
#include <iostream>

Converter::Converter()
{
	_nLimit = 9;
	_matrix = NULL;
	_v = NULL;
	_w = NULL;
	_vPtrs = NULL;
	_matPtrs = NULL;
	_compObject = NULL;
	_comp = NULL;
	_scaleObject = NULL;
	_scale = NULL;
}

void Converter::setupConverter(int count)
{
	_nParam = count;

	double size = sizeof(double) * count * count;
	_matrix = (double *)malloc(size);
	_v = (double *)malloc(size);

	memset(_matrix, '\0', size);
	memset(_v, '\0', size);
	
	_matPtrs = (double **)malloc(sizeof(double **) * count);
	_vPtrs = (double **)malloc(sizeof(double **) * count);
	_w = (double *)malloc(sizeof(double *) * count);
	
	for (int i = 0; i < count; i++)
	{
		_matPtrs[i] = &_matrix[i * count];
		_vPtrs[i] = &_v[i * count];
	}
}

void Converter::setStrategy(RefinementStrategyPtr strategy)
{
	_strategy = strategy;
	setupConverter(_strategy->parameterCount());

	for (int i = 0; i < _strategy->parameterCount(); i++)
	{
		Parameter p = _strategy->getParamObject(i);
		addColumn(p);
	}

	performSVD();
}

void Converter::addColumn(Parameter param)
{
	if (_columns.size() >= _nParam)
	{
		std::cout << "Have used up all column space! Aborting" << std::endl;
		return;
	}

	SVDCol col;
	col.oldParam = param;
	void *o = param.object;
	Getter g = param.getter;

	double default_val = (*g)(o);
	col.start = default_val;
	col.oldParam.start_value = default_val;
	col.inactive = false;
	Param p;
	p.set_value(0);
	col.param = p;
	_columns.push_back(col);
}

void Converter::setCompareFunction(void *obj, CompareParams comp)
{
	_comp = comp;
	_compObject = obj;
}

void Converter::setScaleFunction(void *obj, ScaleParam scale)
{
	_scale = scale;
	_scaleObject = obj;
}

void Converter::scaleColumns()
{
	if (_scale == NULL)
	{
		return;
	}
	std::cout << "Scaling columns..." << std::endl;

	int n = _columns.size();
	for (int j = 0; j < 5; j++)
	{
		for (int i = 0; i < n; i++)
		{
			double scale = (*_scale)(_scaleObject, _columns[i].oldParam);
			if (j == 4)
			{
				std::cout << "Scale " << i << " : " << scale << std::endl;
			}
		}
	}
}

void Converter::compareColumns()
{
	double score = 0;
	int n = _columns.size();
	
//	std::cout << "Comparing columns: " << std::endl;
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				_matPtrs[i][i] = 0.;
				continue;
			}
			
			if (i < j)
			{
				_columns[i].oldParam.step_size *= -1;
			}

			double score = (*_comp)(_compObject, _columns[i].oldParam, 
			                        _columns[j].oldParam);
			
			if (i < j)
			{
				_columns[i].oldParam.step_size *= -1;
			}

			_matPtrs[i][j] = score;
		}
	}
}

void Converter::performSVD()
{
	if (_columns.size() != _nParam)
	{
		std::cout << "Aborting SVD: no match for number of params" << std::endl;
		return;
	}

	scaleColumns();
	
	if (!_compObject)
	{
		return;
	}
	
	compareColumns();
	
	/*
	std::cout << "Pre-SVD results: " << std::endl;
	for (int i = 0; i < _nParam; i++)
	{
		for (int j = 0; j < _nParam; j++)
		{
			std::cout << _matPtrs[i][j] << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	
	size_t dim = _columns.size();
	int success = svdcmp((mat)_matPtrs, dim, dim, (vect)_w, (mat)_vPtrs);
	
	if (!success)
	{
		std::cout << "SVD failure." << std::endl;
		return;
	}

	/*
	std::cout << "Post-SVD results: " << std::endl;
	for (int i = 0; i < _nParam; i++)
	{
		std::cout << _w[i] << ", ";
	}
	std::cout << std::endl << std::endl;
	
	for (int i = 0; i < _nParam; i++)
	{
		for (int j = 0; j < _nParam; j++)
		{
			std::cout << _matPtrs[i][j] << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/

	addParamsToStrategy();
}

double Converter::myScore()
{
	/** convert new parameters into old parameters */
	
	for (int i = 0; i < _columns.size(); i++)
	{
		_columns[i].scratch = _columns[i].start;
	}

	/** run nested evaluation function */
	for (int i = 0; i < _columns.size(); i++)
	{
		for (int j = 0; j < _columns.size(); j++)
		{
			if (_columns[j].inactive)
			{
				continue;
			}

			double val = _columns[j].param.value();
			double step = _columns[j].oldParam.step_size;
			double mod = _matPtrs[i][j];
			step *= mod * val;
			double add = step;

			_columns[i].scratch += add;
		}
	}

	for (int i = 0; i < _columns.size(); i++)
	{
		void *obj = _columns[i].oldParam.object;
		(*_columns[i].oldParam.setter)(obj, _columns[i].scratch);
	}
	
	double score = (*_evalFunc)(_evalObject);

	return score;
}

double Converter::score(void *object)
{
	return static_cast<Converter *>(object)->myScore();
}

void Converter::addParamsToStrategy()
{
	_evalObject = _strategy->getEvaluationObject();
	_evalFunc = _strategy->getEvaluationFunction();
	
	_strategy->setEvaluationFunction(Converter::score, this);
	_strategy->setPartialEvaluation(NULL);
	_strategy->clearParameters();

	double step = 1.0;

	for (int i = 0; i < _nParam; i++)
	{
		if (_w[i] < 0.2)
		{
			_columns[i].inactive = true;
			continue;
		}

		_strategy->addParameter(&_columns[i].param, Param::getValue,
		                       Param::setValue, step, step / 200);
	}
}
