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
#include <vector>
#include "vec3.h"
#include <string.h>
#include "../libccp4/csymlib.h"

struct mat3x3
{
	double vals[9];
};

std::string mat3x3_desc(mat3x3 mat);

mat3x3 mat3x3_inverse(mat3x3 &mat);
mat3x3 mat3x3_from_unit_cell(double a, double b, double c, double alpha, double beta, double gamma);
mat3x3 mat3x3_from_unit_cell(double *unitCell);
void unit_cell_from_mat3x3(mat3x3 mat, double *vals);
mat3x3 make_mat3x3();
mat3x3 mat3x3_from_ccp4(CSym::ccp4_symop symop);

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

void mat3x3_mult_scalar(mat3x3 *mat, double scale);
void mat3x3_scale(mat3x3 *mat, double a, double b, double c);
double mat3x3_length(mat3x3 &mat, int index);
mat3x3 mat3x3_transpose(mat3x3 &mat);
double mat3x3_determinant(mat3x3 &mat);
mat3x3 mat3x3_mult_mat3x3(struct mat3x3 m1, struct mat3x3 m2);
mat3x3 mat3x3_unit_vec_rotation(vec3 axis, double radians);
mat3x3 mat3x3_rotate(double alpha, double beta, double gamma);
mat3x3 mat3x3_ortho_axes(vec3 cVec);
mat3x3 mat3x3_rhbasis(vec3 aVec, vec3 cVec);
mat3x3 mat3x3_closest_rot_mat(vec3 vec1, vec3 vec2, vec3 axis,
                              double *best = NULL);
mat3x3 mat3x3_covariance(std::vector<vec3> points);
mat3x3 mat3x3_make_tensor(mat3x3 &tensify, vec3 &lengths);
mat3x3 mat3x3_subtract_mat3x3(mat3x3 &one, mat3x3 &two);
double mat3x3_abs_sum_all(mat3x3 &mat);

mat3x3 mat3x3_rot_from_angles(double phi, double psi);
mat3x3 mat3x3_from_2d_array(double **values);
void mat3x3_to_2d_array(mat3x3 mat, double ***values);
void free_2d_array(double **values);
double mat3x3_diff_from_identity(mat3x3 &mat, double target = -1);

double mat3x3_rotation_angle(mat3x3 &mat);
vec3 mat3x3_rotation_axis(mat3x3 &mat);

#endif /* defined(__vagabond__mat3x3__) */

