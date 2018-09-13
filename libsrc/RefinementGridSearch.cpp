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
#include <float.h>
#include "CSV.h"
#include <iostream>
#include <iomanip>
#include "FileReader.h"

int RefinementGridSearch::_refine_counter = 0;

void RefinementGridSearch::recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results)
{
	size_t paramCount = parameterCount();
	size_t workingCount = workingList.size();
	Parameter *param = &_params[workingCount];
	double grid_length = param->step_size / param->other_value;

	if (workingCount < paramCount)
	{
		for (int i = -grid_length / 2; i <= (int)(grid_length / 2 + 0.5); i++)
		{
			double mean = referenceList[workingCount];
			/* in refinement grid, step is actually limit... oops */
			double step = param->other_value;
			double value = mean + i * step;

			ParamList extended = workingList;
			extended.push_back(value);
			recursiveEvaluation(referenceList, extended, results);
		}

		return;
	}

	for (int i = 0; i < workingList.size(); i++)
	{
		setValueForParam(i, workingList[i]);
	}

	double result = (*evaluationFunction)(evaluateObject);
	(*results)[workingList] = result;
	reverseResults[result] = workingList;

	orderedParams.push_back(workingList);
	orderedResults.push_back(result);

	reportProgress(result);
}

void RefinementGridSearch::refine()
{
	RefinementStrategy::refine();

	ParamList currentValues;
	CSVPtr csv = CSVPtr(new CSV());

	for (int i = 0; i < parameterCount(); i++)
	{
		double val = getValueForParam(i);
		currentValues.push_back(val);
		csv->addHeader(_params[i].tag);
	}

	csv->addHeader("result");

	recursiveEvaluation(currentValues, ParamList(), &results);

	double minResult = (*evaluationFunction)(evaluateObject);
	ParamList minParams;
	bool changed = false;

	for (ResultMap::iterator it = results.begin(); it != results.end(); it++)
	{
		if (it->second < minResult)
		{
			minResult = it->second;
			minParams = it->first;
			changed = true;
		}

		std::vector<double> result = it->first;
		result.push_back(it->second);

		csv->addEntry(result);
	}

	for (int i = 0; i < minParams.size(); i++)
	{
		double value = currentValues[i];

		if (!_mock && changed)
		{
			value = minParams[i];
		}

		setValueForParam(i, value);
	}

	if (parameterCount() == 2 && _writePNG)
	{
		int stride = _params[0].step_size / _params[0].other_value;

		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = jobName + "_gridsearch_" + i_to_str(_refine_counter);
		plotMap["height"] = "800";
		plotMap["width"] = "800";
		plotMap["xHeader0"] = _params[0].tag;
		plotMap["yHeader0"] = _params[1].tag;
		plotMap["zHeader0"] = "result";

		plotMap["xTitle0"] = _params[0].tag;
		plotMap["yTitle0"] = _params[1].tag;
		plotMap["style0"] = "heatmap";
		plotMap["stride"] = i_to_str(stride);

		csv->plotPNG(plotMap);
	}
	if (_writeCSV)
	{
		csv->writeToFile(jobName + "_gridsearch.csv");
	}

	_refine_counter++;

	finish();
}

