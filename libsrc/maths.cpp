//
//  maths.cpp
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "maths.h"
#include <math.h>
#include <vector>
#include <iostream>

double scale_factor(std::vector<double> &set1, std::vector<double> &set2)
{
	double x_squared = 0;
	double x_y = 0;

	for (int i = 0; i < set1.size(); i++)
	{
		if (set1[i] == set1[i] && set2[i] == set2[i])
		{
			x_squared += set1[i] * set2[i];
			x_y += set2[i] * set2[i];
		}
	}

	double grad = (x_y / x_squared);

	if (grad < 0)
		grad = -1;

	return grad;
}

double r_factor(std::vector<double> &set1, std::vector<double> &set2)
{
	double numerator = 0;
	double denominator = 0;

	for (int i = 0; i < set1.size(); i++)
	{
		if (set1[i] == set1[i] && set2[i] == set2[i])
		{
			double amp_x = (set1[i]);
			double amp_y = (set2[i]);

			numerator += fabs(amp_y - amp_x);
			denominator += fabs((amp_x + amp_y) / 2);
		}
	}

	double rfactor = numerator / denominator;

	return rfactor;
}

double correlation(std::vector<double> &vec1, std::vector<double> &vec2)
{
	double sum_x = 0;
	double sum_y = 0;
	double num = 0;

	if (!vec1.size() || !vec2.size())
	{
		return 0;
	}

	for (int i = 0; i < vec1.size(); i++)
	{
		if (vec1[i] != vec1[i] || vec2[i] != vec2[i])
		{
			continue;
		}

		double addition = vec1[i];
		sum_x += vec1[i];;

		addition = vec2[i];
		sum_y += addition;

		num++;
	}

	double mean_x = sum_x / num;
	double mean_y = sum_y / num;

	if (mean_x != mean_x || mean_y != mean_y)
		return 0;

	double sum_x_y_minus_mean_x_y = 0;
	double sum_x_minus_mean_x_sq = 0;
	double sum_y_minus_mean_y_sq = 0;

	for (int i = 0; i < vec1.size(); i++)
	{
		if (vec1[i] != vec1[i] || vec2[i] != vec2[i])
		{
			continue;
		}

		double addition = (vec1[i] - mean_x) * (vec2[i] - mean_y);
		sum_x_y_minus_mean_x_y += addition;

		addition = pow(vec1[i] - mean_x, 2);
		sum_x_minus_mean_x_sq += addition;

		addition = pow(vec2[i] - mean_y, 2);
		sum_y_minus_mean_y_sq += addition;

	}

	sum_x_y_minus_mean_x_y /= num;
	sum_x_minus_mean_x_sq /= num;
	sum_y_minus_mean_y_sq /= num;

	double r = sum_x_y_minus_mean_x_y
	/ (sqrt(sum_x_minus_mean_x_sq) * sqrt(sum_y_minus_mean_y_sq));
	
	return r;
}

