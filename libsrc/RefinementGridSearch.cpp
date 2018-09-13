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
#include "FileReader.h"

int RefinementGridSearch::_refine_counter = 0;

void RefinementGridSearch::recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results)
{
	size_t paramCount = objects.size();
	size_t workingCount = workingList.size();
	double grid_length = stepSizes[workingCount] / otherValues[workingCount];

	if (workingCount < paramCount)
	{
		if (workingCount == 1)
		{
			//     std::cout << "." << std::flush;
		}

		for (int i = -grid_length / 2; i <= (int)(grid_length / 2 + 0.5); i++)
		{
			double mean = referenceList[workingCount];
			double step = otherValues[workingCount];
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

	for (int i = 0; i < objects.size(); i++)
	{
		Getter getter = getters[i];
		currentValues.push_back((*getter)(objects[i]));
		csv->addHeader(tags[i]);
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
		Setter setter = setters[i];

		if (!_mock && changed)
		{
			(*setter)(objects[i], minParams[i]);
		}
		else
		{
			(*setter)(objects[i], currentValues[i]);
		}
	}

	if (tags.size() == 2 && _writePNG)
	{
		int stride = stepSizes[0] / otherValues[0];

		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = jobName + "_gridsearch_" + i_to_str(_refine_counter);
		plotMap["height"] = "800";
		plotMap["width"] = "800";
		plotMap["xHeader0"] = tags[0];
		plotMap["yHeader0"] = tags[1];
		plotMap["zHeader0"] = "result";

		plotMap["xTitle0"] = tags[0];
		plotMap["yTitle0"] = tags[1];
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

