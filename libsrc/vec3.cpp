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

struct vec3 vec3_subtract_vec3(vec3 &aVec, vec3 &bVec)
{
	struct vec3 vec;
	vec.x = aVec.x - bVec.x;
	vec.y = aVec.y - bVec.y;
	vec.z = aVec.z - bVec.z;

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

void vec3_mult(vec3 *aVec, double mult)
{
	aVec->x *= mult;
	aVec->y *= mult;
	aVec->z *= mult;
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
