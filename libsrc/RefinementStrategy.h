//
//  RefinementStrategy.h
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

#ifndef __cppxfel__RefinementStrategy__
#define __cppxfel__RefinementStrategy__

#include <stdio.h>
#include "shared_ptrs.h"
#include <string>
#include <vector>

typedef enum
{
    MinimizationMethodStepSearch = 0,
    MinimizationMethodNelderMead = 1,
    MinimizationMethodGridSearch = 2,
} MinimizationMethod;


typedef double (*Getter)(void *);
typedef void (*Setter)(void *, double newValue);

class RefinementStrategy
{
protected:
    Getter evaluationFunction;
    Getter finishFunction;
    int maxCycles;
    void *evaluateObject;
    std::string jobName;
    int cycleNum;
    int _changed;
    bool _mock;
    bool _toDegrees;
    bool _silent;

    std::vector<int> couplings;
    std::vector<void *> objects;
    std::vector<Getter> getters;
    std::vector<Setter> setters;
    std::vector<double> stepSizes;
    std::vector<double> otherValues;
    std::vector<std::string> tags;
    std::vector<double> startingValues;
    double startingScore;
    bool _verbose;
    
    void reportProgress(double score);
    void finish();

public:
    RefinementStrategy()
    {
        evaluationFunction = NULL;
        maxCycles = 30;
        cycleNum = 0;
        startingScore = 0;
        _verbose = false;
        _silent = false;
        _changed = -1;
        finishFunction = NULL;
        _mock = false;
        _toDegrees = false;
    };
    
    static RefinementStrategyPtr userChosenStrategy();

    void reportInDegrees()
    {
        _toDegrees = true;
    }
    virtual void refine();
    void resetToInitialParameters();
    
    void addParameter(void *object, Getter getter, Setter setter, double stepSize, double otherValue, std::string tag = "");
    void addCoupledParameter(void *object, Getter getter, Setter setter, double stepSize, double otherValue, std::string tag = "");
    
    void setEvaluationFunction(Getter function, void *evaluatedObject)
    {
        evaluationFunction = function;
        evaluateObject = evaluatedObject;
    }

    void setFinishFunction(Getter finishFunc)
    {
        finishFunction = finishFunc;
    }

    void setVerbose(bool value)
    {
        _verbose = value;
    }

    void setSilent(bool silent)
    {
        _silent = silent;
    }
    
    bool didChange()
    {
        return (_changed == 1);
    }

    void isMock()
    {
        _mock = true;
    }
    
    void setCycles(int num)
    {
        maxCycles = num;
    }
    
    void setJobName(std::string job)
    {
        jobName = job;
    }
    
    void *getEvaluationObject()
    {
        return evaluateObject;
    }
    
    virtual void clearParameters()
    {
        getters.clear();
        setters.clear();
        objects.clear();
        stepSizes.clear();
        otherValues.clear();
        tags.clear();
    }

    int parameterCount()
    {
        return getters.size();
    }
};

#endif /* defined(__cppxfel__RefinementStrategy__) */
