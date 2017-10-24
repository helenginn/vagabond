//
//  RefinementGridSearch.h
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

#ifndef __cppxfel__RefinementGridSearch__
#define __cppxfel__RefinementGridSearch__

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
	ReverseMap reverseResults;
    std::vector<double> orderedResults;
    std::vector<ParamList> orderedParams;
	static int _refine_counter; /* thread care! */

public:
    RefinementGridSearch() : RefinementStrategy()
    {
		gridJumps = 8;
        gridLength = 15;
        cycleNum = 1;
		_writeCSV = false;
    };
    
    void setGridLength(int length)
    {
        gridLength = length;
    }
    
    void setCheckGridNum(int _jumps)
    {
        gridJumps = _jumps;
    }
    
    ResultMap results;
    void recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results);

	std::vector<double> getNextResult(int num)
	{
		ReverseMap::iterator it = reverseResults.begin();

		for (int i = 0; i < num; i++)
		{
			it++;
		}

		return it->second;
	}

	virtual void clearParameters()
    {
        orderedResults.clear();
        orderedParams.clear();
        RefinementStrategy::clearParameters();
    }
    virtual void refine();
};

#endif /* defined(__cppxfel__RefinementGridSearch__) */
