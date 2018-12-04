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

double standard_deviation(std::vector<double> &values) 
{
	double squaredSum = 0;
	double weightSqSum = 0;
	double ave = mean(values);

	for (int i = 0; i < values.size(); i++)
	{
		double value = values[i];

		if (value != value || value == FLT_MAX)
		continue;

		squaredSum += pow(ave - value, 2);

		double weight = 1;

		weightSqSum += weight;
	}

	double stdev = sqrt(squaredSum / weightSqSum);

	return stdev;
}

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

double happiness_coefficient(std::vector<double> xs, std::vector<double> ys)
{
	double intercept, gradient;
	regression_line(xs, ys, &intercept, &gradient);
	
	double dists = 0;
	
	for (int i = 0; i < xs.size(); i++)
	{
		double x = xs[i];
		double y = ys[i];
		
		double closest_x = y * gradient + x - gradient;
		closest_x /= (1 + gradient * gradient);
		double closest_y = y * gradient * gradient + gradient * x + intercept;
		closest_y /= (1 + gradient * gradient);
		
		double dist = (closest_x - x) * (closest_x - x);
		dist += (closest_y - y) * (closest_y - y);
		double expo = exp(-dist);
		
		dists += expo;
	}
	
	dists /= (double)xs.size();
	
	return dists;
}

void regression_line(std::vector<double> xs, std::vector<double> ys,
                     double *intercept, double *gradient)
{
	double sigma_x = 0;
	double sigma_y = 0;
	double sigma_x_y = 0;
	double sigma_x_2 = 0;
	double weight_n = 0;

	for (int i=0; i < xs.size(); i++)
	{
		double x = xs[i];
		double y = ys[i];
		double weight = 1;

		sigma_x += x * weight;
		sigma_y += y * weight;
		sigma_x_y += x * y * weight;
		sigma_x_2 += x * x * weight;
		weight_n += weight;
	}

	double mean_x = sigma_x / weight_n;
	double mean_y = sigma_y / weight_n;

	double sxy = sigma_x_y - sigma_x * sigma_y / weight_n;
	double sxx = sigma_x_2 - pow(sigma_x, 2) / weight_n;

	*gradient = sxy / sxx;
	*intercept = mean_y - *gradient * mean_x;
}

double add_if_gt_zero(std::vector<double> &vec2)
{
	double sum = 0;
	
	for (int i = 0; i < vec2.size(); i++)
	{
		if (vec2[i] < 10e-6)
		{
			continue;
		}
		
		sum++;
	}
	
	return sum;
}

double add_x_if_y(std::vector<double> &vec1, std::vector<double> &vec2,
                  int val)
{
	double sum = 0;
	
	for (int i = 0; i < vec1.size(); i++)
	{
		if (vec2[i] < 1e-2)
		{
			continue;
		}

		if (val == 1 && vec1[i] < 0)
		{
			continue;
		}
		else if (val == -1 && vec1[i] > 0)
		{
			continue;
		}
		
		sum += vec1[i];
		
	}
	
	sum /= vec1.size();
	
	return sum;
}

