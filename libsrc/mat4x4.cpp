//
//  mat4x4.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 02/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "mat4x4.h"
#include <string.h>
#include "vec3.h"
#include "mat3x3.h"
#include <sstream>

mat4x4 make_mat4x4()
{
	mat4x4 mat;
	memset(&mat.vals[0], 0, 16 * sizeof(float));
	mat.vals[0] = 1;
	mat.vals[5] = 1;
	mat.vals[10] = 1;
	mat.vals[15] = 1;

	return mat;
}

/* Auto add 1 at end of vec3 */
vec3 mat4x4_mult_vec(struct mat4x4 mat, struct vec3 vec)
{
	struct vec3 v;

	v.x = mat.vals[0] * vec.x + mat.vals[1] * vec.y + mat.vals[2] * vec.z + mat.vals[3];
	v.y = mat.vals[4] * vec.x + mat.vals[5] * vec.y + mat.vals[6] * vec.z + mat.vals[7];
	v.z = mat.vals[8] * vec.x + mat.vals[9] * vec.y + mat.vals[10] * vec.z + mat.vals[11];

	return v;
}

std::string mat4x4_desc(mat4x4 mat)
{
	std::ostringstream str;
	str << "(" << mat.vals[0] << ", " << mat.vals[1] << ", " << mat.vals[2] << ", " << mat.vals[3] << ",\n";
	str << mat.vals[4] << ", " << mat.vals[5] << ", " << mat.vals[6] << ", " << mat.vals[7] << ",\n";
	str << mat.vals[8] << ", " << mat.vals[9] << ", " << mat.vals[10] << ", " << mat.vals[11] << ",\n";
	str << mat.vals[12] << ", " << mat.vals[13] << ", " << mat.vals[14] << ", " << mat.vals[15] << ")";

	return str.str();
}

mat4x4 mat4x4_frustum(float lr, float bt, float near, float far)
{
	mat4x4 mat = make_mat4x4();

	mat.vals[0] = 1 / lr;
	mat.vals[5] = 1 / bt;
	mat.vals[10] = -2 / (far - near);
	mat.vals[11] = -(far + near) / (far - near);
	mat.vals[15] = 1;

	return mat;
}

mat4x4 mat4x4_mult_mat4x4(mat4x4 l, mat4x4 r)
{
	mat4x4 m;
	m.vals[0] = l.vals[0] * r.vals[0] + l.vals[1] * r.vals[4] + l.vals[2] * r.vals[8] + l.vals[3] * r.vals[12];
	m.vals[1] = l.vals[0] * r.vals[1] + l.vals[1] * r.vals[5] + l.vals[2] * r.vals[9] + l.vals[3] * r.vals[13];
	m.vals[2] = l.vals[0] * r.vals[2] + l.vals[1] * r.vals[6] + l.vals[2] * r.vals[10] + l.vals[3] * r.vals[14];
	m.vals[3] = l.vals[0] * r.vals[3] + l.vals[1] * r.vals[7] + l.vals[2] * r.vals[11] + l.vals[3] * r.vals[15];

	m.vals[4] = l.vals[4] * r.vals[0] + l.vals[5] * r.vals[4] + l.vals[6] * r.vals[8] + l.vals[7] * r.vals[12];
	m.vals[5] = l.vals[4] * r.vals[1] + l.vals[5] * r.vals[5] + l.vals[6] * r.vals[9] + l.vals[7] * r.vals[13];
	m.vals[6] = l.vals[4] * r.vals[2] + l.vals[5] * r.vals[6] + l.vals[6] * r.vals[10] + l.vals[7] * r.vals[14];
	m.vals[7] = l.vals[4] * r.vals[3] + l.vals[5] * r.vals[7] + l.vals[6] * r.vals[11] + l.vals[7] * r.vals[15];

	m.vals[8] = l.vals[8] * r.vals[0] + l.vals[9] * r.vals[4] + l.vals[10] * r.vals[8] + l.vals[11] * r.vals[12];
	m.vals[9] = l.vals[8] * r.vals[1] + l.vals[9] * r.vals[5] + l.vals[10] * r.vals[9] + l.vals[11] * r.vals[13];
	m.vals[10] = l.vals[8] * r.vals[2] + l.vals[9] * r.vals[6] + l.vals[10] * r.vals[10] + l.vals[11] * r.vals[14];
	m.vals[11] = l.vals[8] * r.vals[3] + l.vals[9] * r.vals[7] + l.vals[10] * r.vals[11] + l.vals[11] * r.vals[15];

	m.vals[12] = l.vals[12] * r.vals[0] + l.vals[13] * r.vals[4] + l.vals[14] * r.vals[8] + l.vals[15] * r.vals[12];
	m.vals[13] = l.vals[12] * r.vals[1] + l.vals[13] * r.vals[5] + l.vals[14] * r.vals[9] + l.vals[15] * r.vals[13];
	m.vals[14] = l.vals[12] * r.vals[2] + l.vals[13] * r.vals[6] + l.vals[14] * r.vals[10] + l.vals[15] * r.vals[14];
	m.vals[15] = l.vals[12] * r.vals[3] + l.vals[13] * r.vals[7] + l.vals[14] * r.vals[11] + l.vals[15] * r.vals[15];

	return m;
}

void mat4x4_rotate(mat4x4 *mat, double alpha, double beta, double gamma)
{
	mat3x3 both = mat3x3_rotate(alpha, beta, gamma);
	mat4x4 rot4 = make_mat4x4();

	rot4.vals[0] = both.vals[0];
	rot4.vals[1] = both.vals[1];
	rot4.vals[2] = both.vals[2];

	rot4.vals[4] = both.vals[3];
	rot4.vals[5] = both.vals[4];
	rot4.vals[6] = both.vals[5];

	rot4.vals[8] = both.vals[6];
	rot4.vals[9] = both.vals[7];
	rot4.vals[10] = both.vals[8];

	mat4x4 mult = mat4x4_mult_mat4x4(rot4, *mat);
	*mat = mult;
}

void mat4x4_translate(mat4x4 *mat, vec3 centre)
{
	mat4x4 transMat = make_mat4x4();
	transMat.vals[12] = centre.x;
	transMat.vals[13] = centre.y;
	transMat.vals[14] = centre.z;

	mat4x4 mult = mat4x4_mult_mat4x4(transMat, *mat);
	*mat = mult;
}
