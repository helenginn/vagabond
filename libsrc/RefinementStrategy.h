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
typedef void (*TwoDouble)(void *, double value1, double value2);

typedef struct
{
	void *object;
	Getter getter;
	Setter setter;
	double step_size;
	double other_value;
	double start_value;
	std::string tag;
	int coupled;
} Parameter;

class RefinementStrategy
{
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

	virtual ~RefinementStrategy() {};

	static RefinementStrategyPtr userChosenStrategy();

	void reportInDegrees()
	{
		_toDegrees = true;
	}
	virtual void refine();
	void resetToInitialParameters();

	void addParameter(void *object, Getter getter, Setter setter, 
	                  double stepSize, double otherValue, 
	                  std::string tag = "");
	void addCoupledParameter(void *object, Getter getter, Setter setter, 
	                         double stepSize, double otherValue, 
	                         std::string tag = "");

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
		_params.clear();
	}

	int parameterCount()
	{
		return _params.size();
	}

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

	std::vector<Parameter> _params;
	double startingScore;
	bool _verbose;

	double getValueForParam(int i);
	void setValueForParam(int i, double value);
	void reportProgress(double score);
	void finish();

};

#endif /* defined(__cppxfel__RefinementStrategy__) */
