//
//  maths.cpp
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "vagmaths.h"
#include <cmath>
#include <vector>
#include <iostream>

double weighted_r_factor(std::vector<CoordVal> &vals)
{
	double num = 0;
	double den = 0;
	
	double min_fo = 0;

	for (int i = 0; i < vals.size(); i++)
	{
		double fo = vals[i].fo;
		
		if (min_fo > fo)
		{
			min_fo = fo;
		}
	}
	
	min_fo *= 0.75;

	for (int i = 0; i < vals.size(); i++)
	{
		if (vals[i].weight < 1e-6)
		{
			continue;
		}

		double fo = vals[i].fo - min_fo;
		double fc = vals[i].fc;
		double weight = vals[i].weight;
		
		double add = fabs(fc - fo);
		num += add * weight;
		den += fabs(fo) * weight;
	}
	
	double r = num / den;
	
	return r;
}


double correlation(std::vector<CoordVal> &vals)
{
	double sum_x = 0;
	double sum_y = 0;
	double sum_xx = 0;
	double sum_yy = 0;
	double sum_xy = 0;
	double sum_w = 0;

	for (int i = 0; i < vals.size(); i++)
	{
		if (vals[i].weight < 1e-6)
		{
			continue;
		}

		double weight = vals[i].weight;
		double x = vals[i].fo;
		double y = vals[i].fc;
		
		sum_x += x * weight;
		sum_y += y * weight;
		sum_yy += y * y * weight;
		sum_xx += x * x * weight;
		sum_xy += x * y * weight;
		sum_w += weight;
	}
	
	double top = sum_w * sum_xy - sum_x * sum_y;
	double bottom_left = sum_w * sum_xx - sum_x * sum_x;
	double bottom_right = sum_w * sum_yy - sum_y * sum_y;
	
	double r = top / sqrt(bottom_left * bottom_right);
	
	return r;
}


