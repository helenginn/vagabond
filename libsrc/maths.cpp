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

double scale_factor_by_sum(std::vector<double> &set1, std::vector<double> &set2)
{
	double x_sum = 0;
	double y_sum = 0;

	for (int i = 0; i < set1.size(); i++)
	{
		if (set1[i] == set1[i] && set2[i] == set2[i])
		{
			x_sum += set1[i];
			y_sum += set2[i];
		}
	}

	double grad = (y_sum / x_sum);

	if (grad < 0)
	grad = -1;

	return grad;
}

double scale_factor_cutoff(std::vector<double> &set1, std::vector<double> &set2,
                           double cutoff)
{
	double x_y = 0;
	double y_squared = 0;

	for (int i = 0; i < set1.size(); i++)
	{
		if (set2[i] <= cutoff)
		{
			continue;
		}

		if (set1[i] == set1[i] && set2[i] == set2[i])
		{
			x_y += set1[i] * set2[i];
			y_squared += set2[i] * set2[i];
		}
	}

	double grad = (y_squared / x_y);

	if (grad < 0)
	grad = -1;

	return grad;
}

double scale_factor(std::vector<double> &set1, std::vector<double> &set2)
{
	return scale_factor_cutoff(set1, set2);
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
			denominator += fabs(amp_x);
		}
	}

	double rfactor = numerator / denominator;

	return rfactor;
}

double scaled_r_factor(std::vector<double> &set1, std::vector<double> &set2,
                       double cutoff)
{
	double scale = scale_factor_cutoff(set1, set2, cutoff);

	double numerator = 0;
	double denominator = 0;

	for (int i = 0; i < set1.size(); i++)
	{
		if (set2[i] <= cutoff)
		{
			continue;
		}

		if (set1[i] == set1[i] && set2[i] == set2[i])
		{
			double amp_x = (set1[i]);
			double amp_y = (set2[i]) / scale;

			numerator += fabs(amp_y - amp_x);
			denominator += fabs(amp_x);
		}
	}

	double rfactor = numerator / denominator;

	return rfactor;
}

double mean(std::vector<double> &vec1)
{
	double sum_x = 0;
	double sum_weight = 0;
	for (int i = 0; i < vec1.size(); i++)
	{
		if (vec1[i] != vec1[i])
		{
			continue;
		}

		sum_x += vec1[i];

		sum_weight++;
	}

	return sum_x / sum_weight;
}

double weightedMapScore(std::vector<double> &vec1, std::vector<double> &vec2)
{
	double sum_xy = 0;
	double sum_weight = 0;

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
		double weight = vec2[i];

		sum_xy += addition * weight;
		sum_weight += weight;
	}

	return sum_xy / sum_weight;
}

double correlation(std::vector<double> &vec1, std::vector<double> &vec2,
                   double cutoff)
{
	double sum_x = 0;
	double sum_y = 0;
	double sum_weight = 0;

	if (!vec1.size() || !vec2.size())
	{
		return 0;
	}

	for (int i = 0; i < vec1.size(); i++)
	{
		double weight = 1;

		if (vec1[i] != vec1[i] || vec2[i] != vec2[i] || weight != weight)
		{
			continue;
		}

		if (vec2[i] <= cutoff)
		{
			continue;
		}

		sum_x += vec1[i] * weight;
		sum_y += vec2[i] * weight;

		sum_weight += weight;
	}

	double mean_x = sum_x / sum_weight;
	double mean_y = sum_y / sum_weight;

	if (mean_x != mean_x || mean_y != mean_y)
	return 0;

	double sum_x_y_minus_mean_x_y = 0;
	double sum_x_minus_mean_x_sq = 0;
	double sum_y_minus_mean_y_sq = 0;
	int count = 0;

	for (int i = 0; i < vec1.size(); i++)
	{
		double weight = 1;

		if (vec1[i] != vec1[i] || vec2[i] != vec2[i] || weight != weight)
		{
			continue;
		}

		if (vec2[i] <= cutoff)
		{
			continue;
		}

		double addition = (vec1[i] - mean_x) * (vec2[i] - mean_y);
		sum_x_y_minus_mean_x_y += addition * weight;

		addition = pow(vec1[i] - mean_x, 2);
		sum_x_minus_mean_x_sq += addition * weight;

		addition = pow(vec2[i] - mean_y, 2);
		sum_y_minus_mean_y_sq += addition * weight;
	}

	sum_x_y_minus_mean_x_y /= sum_weight;
	sum_x_minus_mean_x_sq /= sum_weight;
	sum_y_minus_mean_y_sq /= sum_weight;

	double r = sum_x_y_minus_mean_x_y
	/ (sqrt(sum_x_minus_mean_x_sq) * sqrt(sum_y_minus_mean_y_sq));

	return r;
}

/* Produces in real space */
void generateResolutionBins(double minD, double maxD,
                            int binCount, std::vector<double> *bins)
{
	double minRadius = (minD == 0) ? 0 : 1 / minD;
	double maxRadius = 1 / maxD;

	if (maxD <= 0)
	{
		std::cout << "Warning: maximum resolution set to 0. Ignoring.";
		return;
	}

	double maxVolume = pow(maxRadius, 3);
	double minVolume = pow(minRadius, 3);
	double totalVolume = maxVolume - minVolume;

	double eachVolume = totalVolume / binCount;

	double r1 = minRadius;
	double r2 = 0;

	bins->push_back(1 / r1);

	for (int i = 0; i < binCount; i++)
	{
		double r2_cubed = pow(r1, 3) + eachVolume;
		r2 = cbrt(r2_cubed);

		bins->push_back(1 / r2);

		r1 = r2;
	}
}
