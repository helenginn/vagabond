//
//  vec3.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__vec3__
#define __vagabond__vec3__

#include <stdio.h>

struct vec3
{
	double x;
	double y;
	double z;
};

struct vec2
{
	double x;
	double y;
};

double vec3_length(vec3 &vec);
struct vec3 empty_vec3();
inline vec3 make_vec3(double x, double y, double z)
{
	struct vec3 vec;
	vec.x = x;
	vec.y = y;
	vec.z = z;

	return vec;
}

struct vec3 vec3_add_vec3(vec3 &aVec, vec3 &bVec);
struct vec3 vec3_subtract_vec3(vec3 &aVec, vec3 &bVec);
void vec3_mult(vec3 *aVec, double mult);

#endif /* defined(__vagabond__vec3__) */
