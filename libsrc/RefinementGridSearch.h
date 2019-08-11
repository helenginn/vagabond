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

#ifndef __vagabond__RefinementGridSearch__
#define __vagabond__RefinementGridSearch__

#include <stdio.h>
#include <map>
#include "RefinementStrategy.h"

typedef std::vector<double> ParamList;
typedef std::map<ParamList, double> ResultMap;
typedef std::map<double, ParamList> ReverseMap;

class RefinementGridSearch : public RefinementStrategy
{
private:
	int gridLength;
	int gridJumps;
	bool _writeCSV;
	bool _writePNG;
	ReverseMap reverseResults;
	std::vector<double> orderedResults;
	std::vector<ParamList> orderedParams;
	static int _refine_counter; /* thread care! */
	std::vector<double> _array2D;

	double getGridLength(size_t which);

public:
	RefinementGridSearch() : RefinementStrategy()
	{
		gridJumps = 8;
		gridLength = 15;
		cycleNum = 1;
		_writeCSV = false;
		_writePNG = false;
	};

	void setGridLength(int length)
	{
		gridLength = length;
	}

	void setCheckGridNum(int _jumps)
	{
		gridJumps = _jumps;
	}
	
	void setWriteCSV(bool write = true)
	{
		_writeCSV = write;
	}

	void setWritePNG(bool write = true)
	{
		_writePNG = write;
	}
	
	std::vector<double> array2D()
	{
		return _array2D;
	}

	ResultMap results;
	void recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results);

	virtual void clearParameters()
	{
		orderedResults.clear();
		orderedParams.clear();
		RefinementStrategy::clearParameters();
	}
	virtual void refine();
};

#endif /* defined(__vagabond__RefinementGridSearch__) */
