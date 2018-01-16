//
//  Anisotropicator.cpp
//  vagabond
//
//  Created by Helen Ginn on 28/10/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "Anisotropicator.h"
#include "../libica/svdcmp.h"
#include <iostream>

Anisotropicator::Anisotropicator()
{
	_isotropicAverage = 0;
}

void Anisotropicator::setPoints(std::vector<vec3> points)
{
	_points = points;

	_tensor = mat3x3_covariance(points);
	findPrincipleAxes();
}

void Anisotropicator::setTensor(mat3x3 tensor)
{
	_tensor = tensor;
	findPrincipleAxes();
}

double Anisotropicator::anisotropyExtent()
{
	double ave = 0;

	for (int i = 0; i < 3; i++)
	{
		double length = vec3_length(_axes[i]);
		double next = vec3_length(_axes[(i + 1) % 3]);
		ave += fabs(length - next);
	}

	return ave;
}

vec3 Anisotropicator::longestAxis()
{
	int index = 0;
	double longestLength = 0;

	for (int i = 0; i < 3; i++)
	{
		double length = vec3_length(_axes[i]);
		if (length > longestLength)
		{
			longestLength = length;
			index = i;
		}
	}

	return _axes[index];
}

void Anisotropicator::findPrincipleAxes()
{
	mat matrix = NULL;
	double w[3];
	mat3x3_to_2d_array(_tensor, &matrix);
	mat v = NULL;
	mat3x3_to_2d_array(_tensor, &v);

	int success = svdcmp(matrix, 3, 3, (vect)w, (mat)v);

	if (!success)
	{
		return;
	}

	mat3x3 axes = mat3x3_from_2d_array(matrix);
	_axisMatrix = axes;
	axes = mat3x3_transpose(_axisMatrix);
	_isotropicAverage = 0;

	for (int i = 0; i < 3; i++)
	{
		_axes[i] = mat3x3_axis(axes, i);
		vec3_set_length(&_axes[i], sqrt(w[i]));
		_isotropicAverage += w[i] / 3;
	}

	vec3 lengths = {w[0], w[1], w[2]};

	_tensor = mat3x3_make_tensor(_axisMatrix, lengths);
}
