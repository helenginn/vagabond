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
