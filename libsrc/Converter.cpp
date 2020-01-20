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
#include <cstring>
#include <iostream>

Converter::Converter(int count)
{
	double size = sizeof(double) * count * count;
	_matrix = (double *)malloc(size);
	memset(_matrix, '\0', size);

	_nParam = count;
	_compObject = NULL;
	_comp = NULL;
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
	Param p;
	p.set_value(default_val);
	col.param = p;
	_columns.push_back(col);
}

void Converter::setCompareFunction(void *obj, CompareParams comp)
{
	_comp = comp;
	_compObject = obj;
}

void Converter::compareColumns()
{
	double score = 0;
	int n = _columns.size();
	
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (i == j)
			{
				_matrix[i * n + j] = 0;
			}

			double score = (*_comp)(_compObject, _columns[i].oldParam, 
			                        _columns[j].oldParam);
			_matrix[j * n + i] = score;
			_matrix[i * n + j] = -score;
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

	if (_compObject == NULL || _comp == NULL)
	{
		std::cout << "Aborting SVD: comparison objects not set up" << std::endl;
		return;
	}

	compareColumns();
	

}


