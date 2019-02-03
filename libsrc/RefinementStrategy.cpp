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

#include "RefinementGridSearch.h"
#include "RefinementStepSearch.h"
#include "RefinementNelderMead.h"
#include "RefinementStrategy.h"
#include "FileReader.h"
#include <iostream>
#include <iomanip>

RefinementStrategyPtr RefinementStrategy::userChosenStrategy()
{
	int miniMethod = MinimizationMethodNelderMead;
	MinimizationMethod method = (MinimizationMethod)miniMethod;
	RefinementStrategyPtr strategy;

	switch (method) {
		case MinimizationMethodStepSearch:
		strategy = boost::static_pointer_cast<RefinementStrategy>(RefinementStepSearchPtr(new RefinementStepSearch()));
		break;
		case MinimizationMethodNelderMead:
		strategy = boost::static_pointer_cast<RefinementStrategy>(NelderMeadPtr(new RefinementNelderMead()));
		break;
		case MinimizationMethodGridSearch:
		strategy = boost::static_pointer_cast<RefinementStrategy>(RefinementGridSearchPtr(new RefinementGridSearch()));
		break;
		default:
		break;
	}

	int cycles = 20;
	strategy->setCycles(cycles);

	return strategy;
}

void RefinementStrategy::addParameter(void *object, Getter getter, Setter setter, double stepSize, double otherValue, std::string tag, Getter gradient)
{
	Parameter param;
	param.object = object;
	param.getter = getter;
	param.gradient = gradient;
	param.setter = setter;
	param.step_size = stepSize;
	param.other_value = otherValue;

	if (!tag.length())
	{
		tag = "object" + i_to_str((int)parameterCount());
	}

	param.tag = tag;
	param.coupled = 1;

	_params.push_back(param);
}

void RefinementStrategy::addCoupledParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag)
{
	int last = parameterCount() - 1;
	_params[last].coupled++;
	addParameter(object, getter, setter, stepSize, stepConvergence, tag);
	_params[last + 1].coupled++;
}

double RefinementStrategy::getGradientForParam(int i)
{		
	Getter gradient = _params[i].gradient;
	void *object = _params[i].object;
	double grad = (*gradient)(object);
	
	return grad;
}

double RefinementStrategy::getValueForParam(int i)
{
	Getter getter = _params[i].getter;
	void *object = _params[i].object;
	double objectValue = (*getter)(object);

	return objectValue;
}

void RefinementStrategy::setValueForParam(int i, double value)
{
	Setter setter = _params[i].setter;
	void *object = _params[i].object;
	(*setter)(object, value);
}

void RefinementStrategy::refine()
{
	if (!jobName.length())
	{
		jobName = "Refinement procedure for " 
		+ i_to_str((int)parameterCount()) + " objects";
	}

	if (parameterCount() == 0)
	{
		std::cout << " No parameters to refine! Exiting." << std::endl;
		return;
	}

	if (evaluationFunction == NULL || evaluateObject == NULL)
	{
		std::cout << "Please set evaluation function and object." << std::endl;
		return;
	}

	startingScore = (*evaluationFunction)(evaluateObject);
	_prevScore = startingScore;

	for (int i = 0; i < parameterCount(); i++)
	{
		double value = getValueForParam(i);
		_params[i].start_value = value;
	}

	reportProgress(startingScore);
}

void RefinementStrategy::reportProgress(double score)
{
	if (!_verbose || _silent)
	{
		return;
	}

	if (score < _prevScore)
	{
		std::cout << "+" << std::flush;
	}
	else
	{
		std::cout << "=" << std::flush;
	}
	
	cycleNum++;
}

void RefinementStrategy::finish()
{
	if (_verbose)
	{
		std::cout << std::endl;
	}
	double endScore = (*evaluationFunction)(evaluateObject);
	
	if (!parameterCount())
	{
		return;
	}
	
	std::cout << std::setprecision(4);

	if (endScore >= startingScore || endScore != endScore)
	{
		resetToInitialParameters();
		_changed = 0;

		if (!_silent)
		{
			double rad2degscale = (_toDegrees ? rad2deg(1) : 1);
			std::cout << "No change for " << jobName << " ";

			for (int i = 0; i < parameterCount(); i++)
			{
				double value = getValueForParam(i);
				_params[i].changed = 0;
				std::cout << _params[i].tag << "=" << value * rad2degscale <<
				(_toDegrees ? "º" : "") << ", ";
			}

			std::cout << " (" << startingScore << ")" << std::endl;
		}
	}
	else
	{
		double reduction = (startingScore - endScore) / startingScore;
		_improvement = -reduction * 100;

		if (!_silent)
		{std::cout << "Reduction ";
			double rad2degscale = (_toDegrees ? rad2deg(1) : 1);

			if (reduction == reduction)
			{
				std::cout << "by " << std::fixed << 
				-reduction * 100 << "% ";
			}

			std::cout << "for " << jobName << ": ";

			for (int i = 0; i < parameterCount(); i++)
			{
				double value = getValueForParam(i);
				double start = _params[i].start_value;
				
				_params[i].changed = (fabs(start - value) > 1e-4);
				
				std::cout << _params[i].tag << "=" << value * rad2degscale <<
				(_toDegrees ? "°" : "") << ", ";
			}

			std::cout << "(" << startingScore << " to " << 
			endScore << ")" << std::endl;
		}

		_changed = 1;
	}

	cycleNum = 0;

	if (finishFunction != NULL)
	{
		(*finishFunction)(evaluateObject);
	}
}

void RefinementStrategy::resetToInitialParameters()
{
	for (int i = 0; i < parameterCount(); i++)
	{
		double value = _params[i].start_value;
		setValueForParam(i, value);
	}
}
