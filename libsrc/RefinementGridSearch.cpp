//
//  RefinementGridSearch.cpp
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

#include "RefinementGridSearch.h"
#include <float.h>
#include "CSV.h"
#include <iostream>
#include <iomanip>

void RefinementGridSearch::recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results)
{
    size_t paramCount = objects.size();
    size_t workingCount = workingList.size();
    
    if (workingCount < paramCount)
    {
        if (workingCount == 1)
        {
       //     std::cout << "." << std::flush;
        }
        
        for (int i = -gridLength / 2; i <= (int)(gridLength / 2 + 0.5); i++)
        {
            double mean = referenceList[workingCount];
            double step = stepSizes[workingCount];
            double value = mean + i * step;
            
            ParamList extended = workingList;
            extended.push_back(value);
            recursiveEvaluation(referenceList, extended, results);
        }
        
        return;
    }
    
    for (int i = 0; i < workingList.size(); i++)
    {
        Setter setter = setters[i];
        (*setter)(objects[i], workingList[i]);
    }
    
    double result = (*evaluationFunction)(evaluateObject);
    (*results)[workingList] = result;
    
    orderedParams.push_back(workingList);
    orderedResults.push_back(result);
    
    reportProgress(result);
}

void RefinementGridSearch::refine()
{
    RefinementStrategy::refine();
    
    ParamList currentValues;
    CSVPtr csv = CSVPtr(new CSV());

    for (int i = 0; i < objects.size(); i++)
    {
        Getter getter = getters[i];
        currentValues.push_back((*getter)(objects[i]));
        csv->addHeader(tags[i]);
    }
    
    csv->addHeader("result");
    
    recursiveEvaluation(currentValues, ParamList(), &results);

    double minResult = FLT_MAX;
    ParamList minParams;
    
    
    for (ResultMap::iterator it = results.begin(); it != results.end(); it++)
    {
        if (it->second < minResult)
        {
            minResult = it->second;
            minParams = it->first;
        }
        
        std::vector<double> result = it->first;
        result.push_back(it->second);
        
        csv->addEntry(result);
    }
    
    std::cout << "Setting params ";
    
    for (int i = 0; i < minParams.size(); i++)
    {
        Setter setter = setters[i];
        (*setter)(objects[i], minParams[i]);
        
        std::cout << tags[i] << " = " << minParams[i] << ", ";
    }
    
    double val = (*evaluationFunction)(evaluateObject);
    std::cout << "score = " << val << std::endl;
    
    csv->writeToFile(jobName + "_gridsearch.csv");
    
    finish();
}

