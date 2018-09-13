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

#include "RefinementStepSearch.h"
#include "FileReader.h"
#include <float.h>

double RefinementStepSearch::minimizeTwoParameters(int whichParam1, int whichParam2, double *bestScore)
{
	double param_trials1[9];
	double param_trials2[9];
	double param_scores[9];

	Parameter *param1 = &_params[whichParam1];
	Parameter *param2 = &_params[whichParam2];
	
	double *meanStep1 = &param1->step_size;
	double *meanStep2 = &param2->step_size;

	if (*meanStep1 < param1->other_value &&
	    *meanStep2 < param2->other_value)
	{
		return 1;
	}

	int j = 0;
	double param_min_score = *bestScore;
	int param_min_num = 4;

	Getter getter1 = param1->getter;
	Setter setter1 = param1->setter;
	void *object1 = param1->object;

	Getter getter2 = param2->getter;
	Setter setter2 = param2->setter;
	void *object2 = param2->object;

	double bestParam1 = (*getter1)(object1);
	double bestParam2 = (*getter2)(object2);

	for (double i = bestParam1 - *meanStep1; j < 3; i += *meanStep1)
	{
		int l = 0;

		for (double k = bestParam2 - *meanStep2; l < 3; k += *meanStep2)
		{
			(*setter1)(object1, i);
			(*setter2)(object2, k);

			double aScore = (*evaluationFunction)(evaluateObject);

			if (aScore != aScore)
			{
				aScore = FLT_MAX;
			}

			param_scores[j * 3 + l] = aScore;
			param_trials1[j * 3 + l] = i;
			param_trials2[j * 3 + l] = k;

			l++;
		}

		j++;
	}

	for (int i = 0; i < 9; i++)
	if (param_scores[i] < param_min_score)
	{
		param_min_score = param_scores[i];
		param_min_num = i;
	}

	(*setter1)(object1, param_trials1[param_min_num]);
	(*setter2)(object2, param_trials2[param_min_num]);

	if (param_min_num == 4)
	{
		*meanStep1 /= 2;
		*meanStep2 /= 2;
	}

	*bestScore = param_min_score;

	return 0;
}

double RefinementStepSearch::minimizeParameter(int whichParam, double *bestScore)
{
	double param_trials[3];
	double param_scores[3];
	Parameter *param = &_params[whichParam];

	double step = param->step_size;
	if (step < param->other_value)
	{
		return 1;
	}

	Getter getter = param->getter;
	Setter setter = param->setter;
	void *object = param->object;

	int j = 0;
	int param_min_num = 1;

	double bestParam = (*getter)(object);

	if (*bestScore != FLT_MAX)
	{
		param_scores[1] = *bestScore;
	}
	else
	{
		double aScore = (*evaluationFunction)(evaluateObject);
		if (aScore != aScore)
		{
			aScore = FLT_MAX;
		}

		param_scores[1] = aScore;
	}

	param_trials[1] = bestParam;

	for (double i = bestParam - step; j < 3; i += step * 2)
	{
		(*setter)(object, i);

		double aScore = (*evaluationFunction)(evaluateObject);

		if (aScore != aScore)
		{
			aScore = FLT_MAX;
		}

		param_scores[j] = aScore;
		param_trials[j] = i;
		j += 2;
	}

	double param_min_score = param_scores[1];

	for (int i = 0; i < 3; i++)
	if (param_scores[i] < param_min_score)
	{
		param_min_score = param_scores[i];
		param_min_num = i;
	}

	(*setter)(object, param_trials[param_min_num]);

	*bestScore = param_min_score;

	if (param_min_num == 1)
	param->step_size /= 2;

	return 0;
}

void RefinementStepSearch::refine()
{
	RefinementStrategy::refine();

	double bestScore = FLT_MAX;

	for (int i = 0; i < maxCycles; i++)
	{
		bool allFinished = true;

		if (afterCycleObject && afterCycleFunction)
		{
			bestScore = FLT_MAX;
		}

		for (size_t j = 0; j < parameterCount(); j++)
		{
			bool coupled = (_params[j].coupled > 1);

			if (!coupled)
			{
				allFinished *= minimizeParameter(j, &bestScore);
			}
			else
			{
				allFinished *= minimizeTwoParameters(j, j + 1, &bestScore);
				j++;
			}
		}

		if (afterCycleObject && afterCycleFunction)
		{
			(*afterCycleFunction)(afterCycleObject);
		}

		reportProgress(bestScore);

		if (allFinished)
		{
			break;
		}
	}

	finish();
}

