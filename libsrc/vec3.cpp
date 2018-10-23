//
//  vec3.c
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "vec3.h"
#include "string.h"
#include <math.h>
#include <iostream>
#include <sstream>


std::string vec3_desc(vec3 vec)
{
	std::ostringstream str;
	str << "(" << vec.x << ", " << vec.y <<
	", " << vec.z << ")";

	return str.str();
}

struct vec3 empty_vec3()
{
	struct vec3 vec;
	memset(&vec.x, 0, 3 * sizeof(double));

	return vec;
}

double vec3_length(vec3 &vec)
{
	return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

double vec3_sqlength(vec3 &vec)
{
	return (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

void vec3_subtract_from_vec3(vec3 *to, vec3 &from)
{
	to->x -= from.x;
	to->y -= from.y;
	to->z -= from.z;
}

struct vec3 vec3_subtract_vec3(vec3 &to, vec3 &from)
{
	struct vec3 vec;
	vec.x = to.x - from.x;
	vec.y = to.y - from.y;
	vec.z = to.z - from.z;

	return vec;
}

struct vec3 vec3_add_vec3(vec3 &aVec, vec3 &bVec)
{
	struct vec3 vec;
	vec.x = aVec.x + bVec.x;
	vec.y = aVec.y + bVec.y;
	vec.z = aVec.z + bVec.z;

	return vec;
}

double vec3_dot_vec3(vec3 &aVec, vec3 &bVec)
{
	double dot_prod = aVec.x * bVec.x + aVec.y * bVec.y + aVec.z * bVec.z;

	return dot_prod;
}

double vec3_angle_with_vec3(vec3 &aVec, vec3 &bVec)
{
	return acos(vec3_cosine_with_vec3(aVec, bVec));
}

double vec3_cosine_with_vec3(vec3 &aVec, vec3 &bVec)
{
	double dot_prod = aVec.x * bVec.x + aVec.y * bVec.y + aVec.z * bVec.z;
	double vec1_length = vec3_length(aVec);
	double vec2_length = vec3_length(bVec);

	double cosTheta = dot_prod / (vec1_length * vec2_length);

	return cosTheta;
}

vec3 vec3_cross_vec3(vec3 &aVec, vec3 &bVec)
{
	double new_h = aVec.y * bVec.z - aVec.z * bVec.y;
	double new_k = aVec.z * bVec.x - aVec.x * bVec.z;
	double new_l = aVec.x * bVec.y - aVec.y * bVec.x;

	return make_vec3(new_h, new_k, new_l);
}

vec3 make_randomish_axis()
{
	struct vec3 vec;
	vec.x = 2 * (rand() / (double)RAND_MAX) - 1;
	vec.y = 2 * (rand() / (double)RAND_MAX) - 1;
	vec.z = 2 * (rand() / (double)RAND_MAX) - 1;

	vec3_set_length(&vec, 1);

	return vec;
}

double vec3_angle_from_three_points(vec3 &aVec, vec3 &bVec, vec3 &cVec)
{
	vec3 aToB = vec3_subtract_vec3(bVec, aVec);
	vec3 aToC = vec3_subtract_vec3(bVec, cVec);

	return vec3_angle_with_vec3(aToB, aToC);
}

double ewald_wavelength(vec3 &index)
{
	double ewald_radius = index.x * index.x + index.y * index.y
	+ index.z * index.z;

	ewald_radius /= (0 - 2 * index.z);
	double ewald_wavelength = 1 / ewald_radius;

	return ewald_wavelength;
}

bool vec3_near_vec3_box(vec3 &pos1, vec3 &pos2, double tol)
{
	if (fabs(pos1.x - pos2.x) > tol) 
	{
		return false;
	}
	if (fabs(pos1.y - pos2.y) > tol) 
	{
		return false;
	}
	if (fabs(pos1.z - pos2.z) > tol) 
	{
		return false;
	}

	return true;
}


