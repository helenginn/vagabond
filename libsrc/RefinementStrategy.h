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

#ifndef __vagabond__RefinementStrategy__
#define __vagabond__RefinementStrategy__

#include <stdio.h>
#include "shared_ptrs.h"
#include <string>
#include <vector>
#include "Timer.h"

typedef enum
{
	MinimizationMethodStepSearch = 0,
	MinimizationMethodNelderMead = 1,
	MinimizationMethodGridSearch = 2,
} MinimizationMethod;

typedef void (*TwoDouble)(void *, double value1, double value2);

typedef struct
{
	void *object;
	Getter getter;
	Getter gradient;
	Setter setter;
	double step_size;
	double other_value;
	double start_value;
	std::string tag;
	int coupled;
	int changed;
} Parameter;

/** \class RefinementStrategy
 *  \brief Abstract class upon which target function optimisers can be built. 
 **/



class RefinementStrategy
{
public:
	RefinementStrategy()
	{
		_enough = false;
		evaluationFunction = NULL;
		_partial = NULL;
		maxCycles = 30;
		cycleNum = 0;
		startingScore = 0;
		_verbose = false;
		_silent = false;
		_changed = -1;
		finishFunction = NULL;
		_mock = false;
		_improvement = 0;
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
	
	bool changedSignificantly()
	{
		return _enough;
	}

	void addParameter(void *object, Getter getter, Setter setter, 
	                  double stepSize, double otherValue, 
	                  std::string tag = "", Getter gradient = NULL);
	void addCoupledParameter(void *object, Getter getter, Setter setter, 
	                         double stepSize, double otherValue, 
	                         std::string tag = "");

	void setEvaluationFunction(Getter function, void *evaluatedObject)
	{
		evaluationFunction = function;
		evaluateObject = evaluatedObject;
	}

	void setPartialEvaluation(PartialScore function)
	{
		_partial = function;
	}

	void setFinishFunction(Getter finishFunc)
	{
		finishFunction = finishFunc;
	}

	void setVerbose(bool value)
	{
		_verbose = value;
	}

	void setSilent(bool silent = true)
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
	
	Getter getEvaluationFunction()
	{
		return evaluationFunction;
	}

	virtual void clearParameters()
	{
		_changed = false;
		_params.clear();
	}

	int parameterCount()
	{
		return _params.size();
	}
	
	Parameter getParamObject(int i)
	{
		return _params[i];
	}
	
	Parameter *getParamPtr(int i)
	{
		return &_params[i];
	}
	
	void removeParameter(int i)
	{
		_params.erase(_params.begin() + i);
	}
	
	void addParameter(Parameter &param)
	{
		_params.push_back(param);
	}
	
	bool didChange(int i)
	{
		return _params[i].changed;
	}

	double improvement()
	{
		return fabs(_improvement);
	}
	
	void reportResult();
protected:
	double degMult()
	{
		return (_toDegrees ? 180 / M_PI : 1);
	}
	
	Getter evaluationFunction;
	Getter finishFunction;
	PartialScore _partial;
	int maxCycles;
	void *evaluateObject;
	std::string jobName;
	int cycleNum;
	int _changed;
	bool _mock;
	bool _toDegrees;
	bool _silent;
	double _improvement;

	std::vector<Parameter> _params;
	double startingScore;
	double _prevScore;
	bool _verbose;
	bool _enough;

	void findIfSignificant();
	double getGradientForParam(int i);
	double estimateGradientForParam(int i);
	double getValueForParam(int i);
	void setValueForParam(int i, double value);
	void reportProgress(double score);
	void finish();

	Timer _timer;
};

#endif /* defined(__vagabond__RefinementStrategy__) */
