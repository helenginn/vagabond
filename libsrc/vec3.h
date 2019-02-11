//
//  vec3.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__vec3__
#define __vagabond__vec3__

#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <string>

struct vec3
{
	double x;
	double y;
	double z;
};

/** Not really used much here */
struct vec2
{
	double x;
	double y;
};

/** Return the length of the vector by pythagoras */
double vec3_length(vec3 &vec);

/** Returns the square of the vector by pythagoras */
double vec3_sqlength(vec3 &vec);

/** Rescale the length of the vector to a given value by
 * altering internal components */
void vec3_set_length(vec3 *vec, double length);

/** create a fresh vec3 set to (0, 0, 0) */
struct vec3 empty_vec3();

/** create a fresh vec3 with (x, y, z) of your choice */
inline vec3 make_vec3(double x, double y, double z)
{
	struct vec3 vec;
	vec.x = x;
	vec.y = y;
	vec.z = z;

	return vec;
}

/** Not a very serious attempt to just create a randomish
 * unit vector. */
vec3 make_randomish_axis();

/** Create a fresh vec2 set to (x, y) */
inline vec2 make_vec2(double x, double y)
{
	struct vec2 vec;
	vec.x = x;
	vec.y = y;

	return vec;
}

inline bool vec2_less_vec2(vec2 x, vec2 y)
{
	return x.x > y.x;
}

/** Modifies bVec by adding aVec */
inline void vec3_add_to_vec3(vec3 *bVec, vec3 &aVec)
{
	bVec->x += aVec.x;
	bVec->y += aVec.y;
	bVec->z += aVec.z;
}

/** Adds aVec onto bVec and returns a separate vec3 */
struct vec3 vec3_add_vec3(vec3 &aVec, vec3 &bVec);

/** Returns a subtraction vector of the form "to minus from" */
struct vec3 vec3_subtract_vec3(vec3 &to, vec3 &from);

void vec3_subtract_from_vec3(vec3 *to, vec3 &from);

/** Calculates the dot product of two vectors */
double vec3_dot_vec3(vec3 &aVec, vec3 &bVec);

/** Returns the cosine between two vectors */
double vec3_cosine_with_vec3(vec3 &aVec, vec3 &bVec);

/** Returns the angle between two vectors, within the range
 * 0-180 (I think! acosine of the value of the function above.) */
double vec3_angle_with_vec3(vec3 &aVec, vec3 &bVec);

/** Returns the cross product of aVec and bVec */
vec3 vec3_cross_vec3(vec3 &aVec, vec3 &bVec);

/** Calculates the vectors a-to-b and a-to-c and reports the
 * angle between them */
double vec3_angle_from_three_points(vec3 &aVec, vec3 &bVec, vec3 &cVec);

/** If a vector were lying on an Ewald sphere passing through the
 * origin and centred at (0, 0, -1/wavelength), what would the
 * value of wavelength be? */
double ewald_wavelength(vec3 &aVec);

/** Return a string describing the vector in human-readable
 * format */
std::string vec3_desc(vec3 vec);

inline void vec3_min_each(vec3 *minVec, vec3 &aVec)
{
	if (aVec.x < minVec->x) minVec->x = aVec.x;
	if (aVec.y < minVec->y) minVec->y = aVec.y;
	if (aVec.z < minVec->z) minVec->z = aVec.z;
}

inline void vec3_max_each(vec3 *maxVec, vec3 &aVec)
{
	if (aVec.x > maxVec->x) maxVec->x = aVec.x;
	if (aVec.y > maxVec->y) maxVec->y = aVec.y;
	if (aVec.z > maxVec->z) maxVec->z = aVec.z;
}
inline void vec3_mult(vec3 *aVec, double mult)
{
	aVec->x *= mult;
	aVec->y *= mult;
	aVec->z *= mult;
}

inline vec3 vec3_mult(vec3 &vec, double mult)
{
	vec3 aVec;
	aVec.x = mult * aVec.x;
	aVec.y = mult * aVec.y;
	aVec.z = mult * aVec.z;
	
	return aVec;
}


inline void vec3_set_length(vec3 *vec, double length)
{
	double now = vec3_length(*vec);
	vec3_mult(vec, length / now);
}

bool vec3_near_vec3_box(vec3 &pos1, vec3 &pos2, double tol);


#endif /* defined(__vagabond__vec3__) */
