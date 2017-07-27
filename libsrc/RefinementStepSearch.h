//
//  RefinementStepSearch.h
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

#ifndef __cppxfel__RefinementStepSearch__
#define __cppxfel__RefinementStepSearch__

#include <stdio.h>
#include "shared_ptrs.h"
#include "RefinementStrategy.h"

class RefinementStepSearch : public RefinementStrategy
{
private:
    double minimizeParameter(int i, double *bestScore);
    double minimizeTwoParameters(int whichParam1, int whichParam2, double *bestScore);
    
    Getter afterCycleFunction;
    void *afterCycleObject;

public:
    RefinementStepSearch() : RefinementStrategy()
    {
        afterCycleFunction = NULL;
        afterCycleObject = NULL;
    };
    
    void setAfterCycleFunction(Getter function, void *evaluatedObject)
    {
        afterCycleFunction = function;
        afterCycleObject = evaluatedObject;
    }
    
    virtual void refine();
    
};

#endif /* defined(__cppxfel__RefinementStepSearch__) */
