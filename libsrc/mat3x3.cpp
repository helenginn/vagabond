//
//  mat3x3.c
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#define deg2rad(a) ((a)*M_PI/180)
#define rad2deg(a) ((a) / M_PI * 180)

#include "mat3x3.h"
#include <string.h>
#include <sstream>
#include "vec3.h"
#include <math.h>
#include <vector>
#include <iostream>

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
	mat3x3 m = make_mat3x3();

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
	inv.vals[3] = -(d * i - f * g) / det;
	inv.vals[6] = (d * h - e * g) / det;
	inv.vals[1] = - (b * i - c * h) / det;
	inv.vals[4] = (a * i - c * g) / det;
	inv.vals[7] = - (a * h - b * g) / det;
	inv.vals[2] = (b * f - c * e) / det;
	inv.vals[5] = - (a * f - c * d) / det;
	inv.vals[8] = (a * e - b * d) / det;

	return inv;

}

mat3x3 mat3x3_from_unit_cell(double *unitCell)
{
	return mat3x3_from_unit_cell(unitCell[0], unitCell[1],
								 unitCell[2], unitCell[3],
								 unitCell[4], unitCell[5]);
}

void unit_cell_from_mat3x3(mat3x3 mat, double *vals)
{
	vec3 a = mat3x3_axis(mat, 0);
	vec3 b = mat3x3_axis(mat, 1);
	vec3 c = mat3x3_axis(mat, 2);

	double alpha = vec3_angle_with_vec3(b, c);
	double beta = vec3_angle_with_vec3(a, c);
	double gamma = vec3_angle_with_vec3(a, b);

	vals[0] = vec3_length(a);
	vals[1] = vec3_length(b);
	vals[2] = vec3_length(c);
	vals[3] = rad2deg(alpha);
	vals[4] = rad2deg(beta);
	vals[5] = rad2deg(gamma);
}

mat3x3 mat3x3_from_unit_cell(double a, double b, double c, double alpha, double beta, double gamma)
{
	double cosA = cos(deg2rad(alpha));
	double cosB = cos(deg2rad(beta));
	double cosC = cos(deg2rad(gamma));

	if (alpha + beta + gamma > 360)
	{
		std::cout << "Unit cell calculation problem - angles add up to "
		<< alpha + beta + gamma << "º!" << std::endl;
	}

	double sinC = sin(deg2rad(gamma));


	double vol_bit = 1 - cosA * cosA - cosB * cosB - cosC * cosC;
	vol_bit += 2 * cosA * cosB * cosC;
	double volume = a * b * c * sqrt(vol_bit);

	mat3x3 mat;
	mat.vals[0] = a;
	mat.vals[3] = 0;
	mat.vals[6] = 0;
	mat.vals[1] = cosC * b;
	mat.vals[4] = sinC * b;
	mat.vals[7] = 0;
	mat.vals[2] = cosB * c;
	mat.vals[5] = c * (cosA - cosB * cosC) / sinC;
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

mat3x3 mat3x3_rotate(double alpha, double beta, double gamma)
{
	vec3 xAxis = {1, 0, 0};
	vec3 yAxis = {0, 1, 0};
	vec3 zAxis = {0, 0, 1};

	mat3x3 xRot = mat3x3_unit_vec_rotation(xAxis, alpha);
	mat3x3 yRot = mat3x3_unit_vec_rotation(yAxis, beta);
	mat3x3 xyRot = mat3x3_mult_mat3x3(yRot, xRot);
	mat3x3 zRot = mat3x3_unit_vec_rotation(zAxis, gamma);
	mat3x3 xyzRot = mat3x3_mult_mat3x3(zRot, xyRot);

	return xyzRot;
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

	mat.vals[0] = aVec.x;
	mat.vals[1] = cVec.x;
	mat.vals[2] = bVec.x;

	mat.vals[3] = aVec.y;
	mat.vals[4] = cVec.y;
	mat.vals[5] = bVec.y;

	mat.vals[6] = aVec.z;
	mat.vals[7] = cVec.z;
	mat.vals[8] = bVec.z;

	return mat;
}



/* Rotate vector (vec1) around axis (axis) by angle theta. Find value of
 * theta for which the angle between (vec1) and (vec2) is minimised. */
mat3x3 mat3x3_closest_rot_mat(vec3 vec1, vec3 vec2, vec3 axis, double *best)
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

//	double sinSq = tan_theta * tan_theta / (1 + tan_theta * tan_theta);
	double sinSq = pow(sin(theta), 2);
	double cc = sqrt(1 - sinSq);
	double s = sqrt(sinSq);

	/* Now we have two possible solutions, theta or theta+pi
	 * and we need to work out which one. This could potentially be
	 * simplified - do we really need so many cos/sins? maybe check
	 * the 2nd derivative instead? */
	double C = 1 - cc;
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

	if (best != NULL)
	{
		*best = bestAngle;
		return make_mat3x3();
	}
	else
	{
		/* Don't return an identity matrix which has been rotated by
		 * theta around "axis", but do assign it to twizzle. */
		return mat3x3_unit_vec_rotation(axis, bestAngle);
	}
}

vec3 mat3x3_axis(mat3x3 me, int i)
{
	vec3 axis;
	axis = make_vec3(me.vals[0 + i], me.vals[3 + i], me.vals[6 + i]);

	return axis;
}

mat3x3 mat3x3_from_ccp4(CSym::ccp4_symop symop)
{
	mat3x3 matrix = make_mat3x3();
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix.vals[i * 3 + j] = symop.rot[i][j];
		}
	}

	return matrix;
}

std::string mat3x3_desc(mat3x3 mat)
{
	std::ostringstream str;
	str << "(" << mat.vals[0] << ", " << mat.vals[1] << ", " << mat.vals[2] << ",\n";
	str << mat.vals[3] << ", " << mat.vals[4] << ", " << mat.vals[5] << ",\n";
	str << mat.vals[6] << ", " << mat.vals[7] << ", " << mat.vals[8] << ")";

	return str.str();
}

mat3x3 mat3x3_from_2d_array(double **values)
{
	mat3x3 new_mat = make_mat3x3();

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			new_mat.vals[j * 3 + i] = values[j][i];
		}
	}

	return new_mat;
}

void mat3x3_to_2d_array(mat3x3 mat, double ***values)
{
	*values = (double **)malloc(sizeof(double *) * 3);

	for (int j = 0; j < 3; j++)
	{
		(*values)[j] = (double *)malloc(sizeof(double) * 3);

		for (int i = 0; i < 3; i++)
		{
			(*values)[j][i] = mat.vals[j * 3 + i];
		}
	}
}

void free_2d_array(double **values)
{
	for (int i = 0; i < 3; i++)
	{
		free(values[i]);
	}

	free(values);
}

mat3x3 mat3x3_covariance(std::vector<vec3> points)
{
	mat3x3 mat = make_mat3x3();
	memset(mat.vals, 0, sizeof(double) * 9);

	vec3 mean = make_vec3(0, 0, 0);

	for (int i = 0; i < points.size(); i++)
	{
		mean = vec3_add_vec3(mean, points[i]);
	}

	vec3_mult(&mean, 1 / (double)points.size());

	for (int i = 0; i < points.size(); i++)
	{
		points[i] = vec3_subtract_vec3(points[i], mean);
	}

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int k = 0; k < points.size(); k++)
			{
				double add = *(&points[k].x + i) * *(&points[k].x + j);
				mat.vals[j * 3 + i] += add;
			}
		}
	}

	for (int i = 0; i < 9; i++)
	{
		mat.vals[i] /= (double)points.size();
	}

	return mat;
}

mat3x3 mat3x3_rot_from_angles(double phi, double psi)
{
	vec3 zAxis = {0, 0, 1};
	vec3 xAxis = {1, 0, 0};
	/* we shall rotate phi round the y axis for starters */
	mat3x3 mat = mat3x3_rotate(0, phi, 0);
	vec3 rotAxis = mat3x3_mult_vec(mat, xAxis);

	/* we need to remake our next rotation axis as the cross product */
	vec3 newAxis = vec3_cross_vec3(rotAxis, zAxis);

	mat3x3 secondRot = mat3x3_unit_vec_rotation(newAxis, psi);

	return mat3x3_mult_mat3x3(secondRot, mat);
}

mat3x3 mat3x3_make_tensor(mat3x3 &tensify, vec3 &lengths)
{
	mat3x3 transpose = mat3x3_inverse(tensify);
	mat3x3 scaling = make_mat3x3();
	scaling.vals[0] = lengths.x;
	scaling.vals[4] = lengths.y;
	scaling.vals[8] = lengths.z;
	mat3x3 combo1 = mat3x3_mult_mat3x3(scaling, transpose);
	mat3x3 combo2 = mat3x3_mult_mat3x3(tensify, combo1);
	
	return combo2;
}

double mat3x3_diff_from_identity(mat3x3 &mat, double target)
{
	double diff = 0;
	if (target < 0)
	{
		target = (mat.vals[0] + mat.vals[4] + mat.vals[8]) / 3;
	}

//	std::cout << "---" << std::endl;

	for (int i = 0; i < 9; i++)
	{
		double this_target = (i % 4 == 0) ? target : 0;
		double add = fabs(mat.vals[i] - this_target);
	//	std::cout << this_target << " " << mat.vals[i] << std::endl;

		diff += add;
	}

	return diff;
}

void mat3x3_mult_scalar(mat3x3 *mat, double scale)
{
	for (int i = 0; i < 9; i++)
	{
		mat->vals[i] *= scale;
	}
}

double mat3x3_rotation_angle(mat3x3 &mat)
{
    double angle = 0;
    double a = mat.vals[0];
    double b = mat.vals[5];
    double c = mat.vals[8];

    double cosTheta = (a + b + c - 1) / 2;
    double theta = acos(cosTheta);

    return theta;
}

vec3 mat3x3_rotation_axis(mat3x3 &mat)
{
    double angle = mat3x3_rotation_angle(mat);
    double cosangle = cos(angle);
    double x = sqrt((mat.vals[0] - cosangle) / (1 - cosangle));
    double y = sqrt((mat.vals[5] - cosangle) / (1 - cosangle));
    double z = sqrt((mat.vals[8] - cosangle) / (1 - cosangle));

    return make_vec3(x, y, z);
}
