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

/** Create a C++ string describing a mat3x3 matrix in
 * human readable format. */
std::string mat3x3_desc(mat3x3 mat);

/** Simple matrix inversion */
mat3x3 mat3x3_inverse(mat3x3 &mat);

/** Calculate basis vectors for a unit cell and organise them
 * into a matrix, using six individual variables */
mat3x3 mat3x3_from_unit_cell(double a, double b, double c, double alpha, double beta, double gamma);

/** Calculate basis vectors for a unit cell and organise them
 * into a matrix, using a pointer to a double[6] array */
mat3x3 mat3x3_from_unit_cell(double *unitCell);

/** Calculate a unit cell from the basis vectors in a matrix */
void unit_cell_from_mat3x3(mat3x3 mat, double *vals);

/** Return an identity matrix */
mat3x3 make_mat3x3();

/** mat3x3 matrix from CCP4 symmetry operation as a 
 * rotation matrix */
mat3x3 mat3x3_from_ccp4(CSym::ccp4_symop symop);

/** Inline method (mildly faster) to multiply a vector by
 * a matrix, by supplying pointer to vector to change */
inline void mat3x3_mult_vec(struct mat3x3 mat, struct vec3 *vec)
{
	struct vec2 v;

	v.x = mat.vals[0] * vec->x + mat.vals[1] * vec->y + mat.vals[2] * vec->z;
	v.y = mat.vals[3] * vec->x + mat.vals[4] * vec->y + mat.vals[5] * vec->z;
	vec->z = mat.vals[6] * vec->x + mat.vals[7] * vec->y + mat.vals[8] * vec->z;

	memcpy(vec, &v.x, sizeof(double) * 2);
}

/** Grab the i'th basis vector from a matrix (column 0, 1, or 2) */
vec3 mat3x3_axis(mat3x3 &me, int i);

/** Add contents of const into alter */
void mat3x3_add_mat3x3(mat3x3 *alter, mat3x3 &_const);

/** Matrix multiplied by vector and returned as a fresh vector.
 * Consider using the inline function above if you don't need to
 * keep the old one. */
vec3 mat3x3_mult_vec(struct mat3x3 mat, struct vec3 vec);

/** Multiply every basis vector in the matrix by a scalar value */
void mat3x3_mult_scalar(mat3x3 *mat, double scale);

/** Multiply each basis vector by its own scalar value. */
void mat3x3_scale(mat3x3 *mat, double a, double b, double c);

/** Calculate the length of the i'th basis vector in a matrix */
double mat3x3_length(mat3x3 &mat, int index);

/** Transpose the matrix and return as a new one, so swapping
 * rows and columns. Also quick way of inverting a rotation
 * matrix */
mat3x3 mat3x3_transpose(mat3x3 &mat);

/** Calculate the determinant of the matrix */
double mat3x3_determinant(mat3x3 &mat);

/** Calculate the trace of the matrix */
double mat3x3_trace(mat3x3 &mat);

/** Volume of unit cell described by matrix mat */
double mat3x3_volume(mat3x3 mat);

/** Multiply two matrices and return the result */
mat3x3 mat3x3_mult_mat3x3(struct mat3x3 m1, struct mat3x3 m2);

/** Calculate a right-handed rotation matrix which is a rotation 
 * around an axis by a given degree in radians. */
mat3x3 mat3x3_unit_vec_rotation(vec3 axis, double radians);

/** Calculate a right-handed rotation matrix specified by
 * rotations around each axis x, y, z - mostly suitable for small
 * angles when the small-angle approximation holds. */
mat3x3 mat3x3_rotate(double alpha, double beta, double gamma);

/** Quickly finds any orthonormal matrix where the last basis
 * vector is specified by cVec, and aVec and bVec complete some
 * kind of axis. Handedness not guaranteed. Everything comes back
 * with unit vector bases. */
mat3x3 mat3x3_ortho_axes(vec3 cVec);

/** Calculates a right-handed matrix based on two, not necessarily
 * right-angled basis vectors. Everything comes back with unit
 * vector bases */
mat3x3 mat3x3_rhbasis(vec3 aVec, vec3 cVec);

/** Find the rotation matrix which maps vec1 as closely as possible
 * onto vec2 by rotating around 'axis'. "best" returns the final
 * closest angle between vec2 and vec1 after rotation. If unity
 * is set to true, then you must guarantee that the incoming
 * vectors are unit vectors, in which case those calculations
 * are skipped and it's a bit quicker. */
mat3x3 mat3x3_closest_rot_mat(vec3 vec1, vec3 vec2, vec3 axis,
                              double *best = NULL, bool unity = false);

/** Supply a C++ std::vector of vec3 positions and this will
 *  return the covariance matrix of this list. */
mat3x3 mat3x3_covariance(std::vector<vec3> points);

/** Multiplies each basis vector of mat3x3 &tensify by the appropriate
 * length in vec3 &lengths, then multiplies by its transpose... might
 * be buggy */
mat3x3 mat3x3_make_tensor(mat3x3 &tensify, vec3 &lengths);

/* Set each axis down columns to length of 1 */
void mat3x3_vectors_to_unity(mat3x3 *mat);

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

