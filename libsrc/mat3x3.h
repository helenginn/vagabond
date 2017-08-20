//
//  mat3x3.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__mat3x3__
#define __vagabond__mat3x3__

#include <stdio.h>
#include "vec3.h"
#include <string.h>

struct mat3x3
{
	double vals[9];
};


mat3x3 mat3x3_inverse(mat3x3 &mat);
mat3x3 mat3x3_from_unit_cell(double a, double b, double c, double alpha, double beta, double gamma);
mat3x3 make_mat3x3();

inline void mat3x3_mult_vec(struct mat3x3 mat, struct vec3 *vec)
{
	struct vec2 v;

	v.x = mat.vals[0] * vec->x + mat.vals[1] * vec->y + mat.vals[2] * vec->z;
	v.y = mat.vals[3] * vec->x + mat.vals[4] * vec->y + mat.vals[5] * vec->z;
	vec->z = mat.vals[6] * vec->x + mat.vals[7] * vec->y + mat.vals[8] * vec->z;

	memcpy(vec, &v.x, sizeof(double) * 2);
}

vec3 mat3x3_axis(mat3x3 me, int i);

vec3 mat3x3_mult_vec(struct mat3x3 mat, struct vec3 vec);

void mat3x3_scale(mat3x3 *mat, double a, double b, double c);
double mat3x3_length(mat3x3 &mat, int index);
mat3x3 mat3x3_transpose(mat3x3 &mat);
double mat3x3_determinant(mat3x3 &mat);
mat3x3 mat3x3_mult_mat3x3(struct mat3x3 m1, struct mat3x3 m2);
mat3x3 mat3x3_unit_vec_rotation(vec3 axis, double radians);
mat3x3 mat3x3_ortho_axes(vec3 cVec);
mat3x3 mat3x3_rhbasis(vec3 aVec, vec3 bVec);
mat3x3 mat3x3_closest_rot_mat(vec3 vec1, vec3 vec2, vec3 axis,
							  double *best = NULL);

#endif /* defined(__vagabond__mat3x3__) */
