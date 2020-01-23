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
	double start;
	double scratch;
} SVDCol;

typedef double (*CompareParams)(void *obj, Parameter &p1, Parameter &p2);
typedef double (*ScaleParam)(void *obj, Parameter &p1);

class Converter
{
public:
	Converter();
	
	void setCompareFunction(void *obj, CompareParams comp);
	void setScaleFunction(void *obj, ScaleParam comp);
	void setStrategy(RefinementStrategyPtr strategy);

	static double score(void *object);
private:
	void setupConverter(int count);
	void addParamsToStrategy();
	void addColumn(Parameter param);
	double myScore();
	void scaleColumns();
	void compareColumns();
	void performSVD();
	/** matrix will contain param-to-param correlations */
	double *_matrix;
	double *_v;

	double **_matPtrs;
	double **_vPtrs;
	double *_w;

	std::vector<SVDCol> _columns;

	/* we copy the evaluation object/function from the strategy and
	 * then act as an interface, converting parameters via SVD
	 * calculations*/

	void *_evalObject;
	Getter _evalFunc;

	RefinementStrategyPtr _strategy;
	void *_compObject;
	CompareParams _comp;

	void *_scaleObject;
	ScaleParam _scale;

	int _nParam;
	int _nLimit;
};

#endif
