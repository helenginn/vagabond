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

double scale_factor(std::vector<double> &set1, std::vector<double> &set2);
double r_factor(std::vector<double> &set1, std::vector<double> &set2);
double correlation(std::vector<double> &vec1, std::vector<double> &vec2);

typedef double (*two_dataset_op)(std::vector<double>&, std::vector<double>&);

inline double normal_distribution(double x, double sigma)
{
	double power = 0 - pow((x), 2) / (2 * sigma * sigma);
	double exp = pow(M_E, power);

	double denominator = sigma * sqrt(2 * M_PI);

	return exp / denominator;
}

#endif /* defined(__vagabond__maths__) */
