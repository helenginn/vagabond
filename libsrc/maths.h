//
//  maths.h
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__maths__
#define __vagabond__maths__

#include <stdio.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <float.h>

double scale_factor_cutoff(std::vector<double> &set1, std::vector<double> &set2,
                           double cutoff = -FLT_MAX);
double scale_factor(std::vector<double> &set1, std::vector<double> &set2);
double scale_factor_by_sum(std::vector<double> &set1, std::vector<double> &set2);

double scaled_r_factor(std::vector<double> &set1, std::vector<double> &set2,
                       double cutoff = -FLT_MAX);
double r_factor(std::vector<double> &set1, std::vector<double> &set2);
double weightedMapScore(std::vector<double> &set1, std::vector<double> &set2);
double correlation(std::vector<double> &vec1, std::vector<double> &vec2,
                   double cutoff = -FLT_MAX);
double happiness_coefficient(std::vector<double> xs, std::vector<double> ys);
double mean(std::vector<double> &vec1);
/* second vec2 is ignored */
inline double two_dataset_mean(std::vector<double> &vec1, std::vector<double> &vec2)
{
	return mean(vec1);
}
double add_x_if_y(std::vector<double> &vec1, std::vector<double> &vec2,
                  int val);
double add_if_gt_zero(std::vector<double> &vec2);
double standard_deviation(std::vector<double> &values);
void regression_line(std::vector<double> xs, std::vector<double> ys,
                     double *intercept, double *gradient);

typedef double (*two_dataset_op)(std::vector<double>&, std::vector<double>&);

/* Dstar */
void generateResolutionBins(double minD, double maxD,
                            int binCount, std::vector<double> *bins);

inline double normal_distribution(double x, double sigma)
{
	double power = 0 - pow((x), 2) / (2 * sigma * sigma);
	double exp = pow(M_E, power);

	double denominator = sigma * sqrt(2 * M_PI);

	return exp / denominator;
}

inline double random_norm_dist(double x, double sigma)
{
	const int tries = 10;
	double total = 0;

	for (int i = 0; i < tries; i++)
	{
		total += rand() / (double)RAND_MAX;
	}

	total -= tries * 0.5;

	total *= sigma;
	total += x;

	return total;
}

#endif /* defined(__vagabond__maths__) */
