//
//  maths.cpp
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "maths.h"
#include <math.h>

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