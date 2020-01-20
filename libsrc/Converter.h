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

#ifndef __vagabond__converter__
#define __vagabond__converter__

#include "Param.h"
#include "RefinementStrategy.h"

/** \brief Allows conversion between parameter sets */

typedef struct
{
	Parameter oldParam; /** old column information */
	Param param;    	/** new parameter to assign SVD version */
} SVDCol;

typedef double (*CompareParams)(void *obj, Parameter &p1, Parameter &p2);

class Converter
{
public:
	Converter(int count);
	
	void addColumn(Parameter param);
	void setCompareFunction(void *obj, CompareParams comp);

	void performSVD();
private:
	void compareColumns();
	/** matrix will contain param-to-param correlations */
	double *_matrix;
	std::vector<SVDCol> _columns;

	void *_compObject;
	CompareParams _comp;
	int _nParam;
};

#endif
