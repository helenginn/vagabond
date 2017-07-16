//
//  mat3x3.c
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#define deg2rad(a) ((a)*M_PI/180)

#include "mat3x3.h"
#include <string.h>
#include "vec3.h"
#include <math.h>

struct mat3x3 make_mat3x3()
{
	struct mat3x3 mat;
	memset(mat.vals, 0, 9 * sizeof(double));

	mat.vals[0] = 1;
	mat.vals[4] = 1;
	mat.vals[8] = 1;
	
	return mat;
}

void mat3x3_mult_vec(struct mat3x3 mat, struct vec3 *vec)
{
	struct vec3 v;

	v.x += mat.vals[0] * vec->x + mat.vals[1] * vec->y + mat.vals[2] * vec->z;
	v.y += mat.vals[3] * vec->x + mat.vals[4] * vec->y + mat.vals[5] * vec->z;
	v.z += mat.vals[6] * vec->x + mat.vals[7] * vec->y + mat.vals[8] * vec->z;

	memcpy(vec, &v.x, sizeof(double) * 3);
}

double mat3x3_determinant(mat3x3 &mat)
{
	double a = mat.vals[0];
	double b = mat.vals[1];
	double c = mat.vals[2];
	double d = mat.vals[3];
	double e = mat.vals[4];
	double f = mat.vals[5];
	double g = mat.vals[6];
	double h = mat.vals[7];
	double i = mat.vals[8];

	double det = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;

	return det;
}


mat3x3 mat3x3_inverse(mat3x3 &mat)
{
	double a = mat.vals[0];
	double b = mat.vals[1];
	double c = mat.vals[2];
	double d = mat.vals[3];
	double e = mat.vals[4];
	double f = mat.vals[5];
	double g = mat.vals[6];
	double h = mat.vals[7];
	double i = mat.vals[8];

	double det = mat3x3_determinant(mat);
	mat3x3 inv;

	inv.vals[0] = (e * i - f * h) / det;
	inv.vals[1] = -(d * i - f * g) / det;
	inv.vals[2] = (d * h - e * g) / det;
	inv.vals[3] = - (b * i - c * h) / det;
	inv.vals[4] = (a * i - c * g) / det;
	inv.vals[5] = - (a * h - b * g) / det;
	inv.vals[6] = (b * f - c * e) / det;
	inv.vals[7] = - (a * f - c * d) / det;
	inv.vals[8] = (a * e - b * d) / det;

	return inv;

}

mat3x3 mat3x3_from_unit_cell(double a, double b, double c, double alpha, double beta, double gamma)
{
	double cosA = cos(deg2rad(alpha));
	double cosB = cos(deg2rad(beta));
	double cosC = cos(deg2rad(gamma));

	double sinC = sin(deg2rad(gamma));

	double volume = a * b * c * sqrt(1 - cosA * cosA - cosB * cosB - cosC * cosC + 2 * cosA * cosB * cosC);

	mat3x3 mat;
	mat.vals[0] = a;
	mat.vals[1] = 0;
	mat.vals[2] = 0;
	mat.vals[3] = cosC * b;
	mat.vals[4] = sinC * b;
	mat.vals[5] = 0;
	mat.vals[6] = cosB * c;
	mat.vals[7] = c * (cosA - cosB * cosC) / sinC;
	mat.vals[8] = volume / (a * b * sinC);

	return mat;
}

mat3x3 mat3x3_transpose(mat3x3 &mat)
{
	mat3x3 new_mat = make_mat3x3();
	new_mat.vals[0] = mat.vals[0];
	new_mat.vals[1] = mat.vals[3];
	new_mat.vals[2] = mat.vals[6];
	new_mat.vals[3] = mat.vals[1];
	new_mat.vals[4] = mat.vals[4];
	new_mat.vals[5] = mat.vals[7];
	new_mat.vals[6] = mat.vals[2];
	new_mat.vals[7] = mat.vals[5];
	new_mat.vals[8] = mat.vals[8];

	return new_mat;
}

void mat3x3_scale(mat3x3 *inverse, double a, double b, double c)
{
	inverse->vals[0] *= a;
	inverse->vals[1] *= a;
	inverse->vals[2] *= a;

	inverse->vals[3] *= b;
	inverse->vals[4] *= b;
	inverse->vals[5] *= b;

	inverse->vals[6] *= c;
	inverse->vals[7] *= c;
	inverse->vals[8] *= c;
}

double mat3x3_length(mat3x3 &mat, int index)
{
	double sqLength = mat.vals[index * 3 + 0] * mat.vals[index * 3 + 0]
	+ mat.vals[index * 3 + 1] * mat.vals[index * 3 + 1]
	+ mat.vals[index * 3 + 2] * mat.vals[index * 3 + 2];

	return sqrt(sqLength);
}
