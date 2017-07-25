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

vec3 mat3x3_mult_vec(struct mat3x3 mat, struct vec3 vec)
{
	struct vec3 v;

	v.x = mat.vals[0] * vec.x + mat.vals[1] * vec.y + mat.vals[2] * vec.z;
	v.y = mat.vals[3] * vec.x + mat.vals[4] * vec.y + mat.vals[5] * vec.z;
	v.z = mat.vals[6] * vec.x + mat.vals[7] * vec.y + mat.vals[8] * vec.z;

	return v;
}

mat3x3 mat3x3_mult_mat3x3(struct mat3x3 m1, struct mat3x3 m2)
{
	struct mat3x3 m;

	m.vals[0] = m1.vals[0] * m2.vals[0] + m1.vals[1] * m2.vals[3] + m1.vals[2] * m2.vals[6];
	m.vals[1] = m1.vals[0] * m2.vals[1] + m1.vals[1] * m2.vals[4] + m1.vals[2] * m2.vals[7];
	m.vals[2] = m1.vals[0] * m2.vals[2] + m1.vals[1] * m2.vals[5] + m1.vals[2] * m2.vals[8];

	m.vals[3] = m1.vals[3] * m2.vals[0] + m1.vals[4] * m2.vals[3] + m1.vals[5] * m2.vals[6];
	m.vals[4] = m1.vals[3] * m2.vals[1] + m1.vals[4] * m2.vals[4] + m1.vals[5] * m2.vals[7];
	m.vals[5] = m1.vals[3] * m2.vals[2] + m1.vals[4] * m2.vals[5] + m1.vals[5] * m2.vals[8];

	m.vals[6] = m1.vals[6] * m2.vals[0] + m1.vals[7] * m2.vals[3] + m1.vals[8] * m2.vals[6];
	m.vals[7] = m1.vals[6] * m2.vals[1] + m1.vals[7] * m2.vals[4] + m1.vals[8] * m2.vals[7];
	m.vals[8] = m1.vals[6] * m2.vals[2] + m1.vals[7] * m2.vals[5] + m1.vals[8] * m2.vals[8];

	return m;
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


mat3x3 mat3x3_unit_vec_rotation(vec3 axis, double radians)
{
	mat3x3 mat;

	double x = axis.x;
	double x2 = axis.x * axis.x;

	double y = axis.y;
	double y2 = axis.y * axis.y;

	double z = axis.z;
	double z2 = axis.z * axis.z;

	double cosa = cos(radians);
	double sina = sin(radians);

	mat.vals[0] = cosa + x2 * (1 - cosa);
	mat.vals[1] = x * y * (1 - cosa) - z * sina;
	mat.vals[2] = x * z * (1 - cosa) + y * sina;

	mat.vals[3] = y * x * (1 - cosa) + z * sina;
	mat.vals[4] = cosa + y2 * (1 - cosa);
	mat.vals[5] = z * y * (1 - cosa) - x * sina;

	mat.vals[6] = z * x * (1 - cosa) - y * sina;
	mat.vals[7] = z * y * (1 - cosa) + x * sina;
	mat.vals[8] = cosa + z2 * (1 - cosa);

	return mat;
}


mat3x3 mat3x3_ortho_axes(vec3 cVec)
{
	vec3_set_length(&cVec, 1);
	vec3 bVec = make_vec3(1, 0, cVec.x / cVec.z);
	vec3_set_length(&bVec, 1);
	vec3 aVec = vec3_cross_vec3(bVec, cVec);

	mat3x3 mat;
	memcpy(&mat.vals[0], &aVec, 3 * sizeof(double));
	memcpy(&mat.vals[3], &bVec, 3 * sizeof(double));
	memcpy(&mat.vals[6], &cVec, 3 * sizeof(double));

	return mat;
}

mat3x3 mat3x3_rhbasis(vec3 aVec, vec3 bVec)
{
	vec3_set_length(&aVec, 1);
	vec3_set_length(&bVec, 1);
	vec3 cVec = vec3_cross_vec3(aVec, bVec);
	vec3_set_length(&cVec, 1);

	mat3x3 mat;
	memcpy(&mat.vals[0], &aVec, 3 * sizeof(double));
	memcpy(&mat.vals[3], &cVec, 3 * sizeof(double));
	memcpy(&mat.vals[6], &bVec, 3 * sizeof(double));

	return mat3x3_transpose(mat);
}



/* Rotate vector (vec1) around axis (axis) by angle theta. Find value of
 * theta for which the angle between (vec1) and (vec2) is minimised. */
mat3x3 mat3x3_closest_rot_mat(vec3 vec1, vec3 vec2, vec3 axis)
{
	/* Let's have unit vectors */
	vec3_set_length(&vec1, 1);
	vec3_set_length(&vec2, 1);
	vec3_set_length(&axis, 1);

	/* Redeclaring these to try and maintain readability and
	 * check-ability against the maths I wrote down */
	double a = vec2.x; double b = vec2.y; double c = vec2.z;
	double p = vec1.x; double q = vec1.y; double r = vec1.z;
	double x = axis.x; double y = axis.y; double z = axis.z;

	/* Components in handwritten maths online when I upload it */
	double A = a*(p*x*x - p + x*y*q + x*z*r) +
	b*(p*x*y + q*y*y - q + r*y*z) +
	c*(p*x*z + q*y*z + r*z*z - r);

	double B = a*(y*r - z*q) + b*(p*z - r*x) + c*(q*x - p*y);

	double tan_theta = - B / A;
	double theta = atan(tan_theta);

	/* Now we have two possible solutions, theta or theta+pi
	 * and we need to work out which one. This could potentially be
	 * simplified - do we really need so many cos/sins? maybe check
	 * the 2nd derivative instead? */
	double cc = cos(theta);
	double C = 1 - cc;
	double s = sin(theta);
	double occ = -cc;
	double oC = 1 - occ;
	double os = -s;

	double pPrime = (x*x*C+cc)*p + (x*y*C-z*s)*q + (x*z*C+y*s)*r;
	double qPrime = (y*x*C+z*s)*p + (y*y*C+cc)*q + (y*z*C-x*s)*r;
	double rPrime = (z*x*C-y*s)*p + (z*y*C+x*s)*q + (z*z*C+cc)*r;

	double pDbPrime = (x*x*oC+occ)*p + (x*y*oC-z*os)*q + (x*z*oC+y*os)*r;
	double qDbPrime = (y*x*oC+z*os)*p + (y*y*oC+occ)*q + (y*z*oC-x*os)*r;
	double rDbPrime = (z*x*oC-y*os)*p + (z*y*oC+x*os)*q + (z*z*oC+occ)*r;

	double cosAlpha = pPrime * a + qPrime * b + rPrime * c;
	double cosAlphaOther = pDbPrime * a + qDbPrime * b + rDbPrime * c;

	int addPi = (cosAlphaOther > cosAlpha);
	double bestAngle = theta + addPi * M_PI;

	/* Don't return an identity matrix which has been rotated by
	 * theta around "axis", but do assign it to twizzle. */
	return mat3x3_unit_vec_rotation(axis, bestAngle);
}

vec3 axis(mat3x3 me, int i)
{
	vec3 axis;
	memcpy(&axis, &me.vals[i * 3], 3 * sizeof(double));

	return axis;
}
