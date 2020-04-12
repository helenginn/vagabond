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

void mat4x4_mult_scalar(struct mat4x4 *mat, double scalar)
{
	for (int i = 0; i < 16; i += 4)
	{
		mat->vals[i] *= scalar;
		mat->vals[i+1] *= scalar;
		mat->vals[i+2] *= scalar;
	}
}

/* Auto add 1 at end of vec3 */
vec3 mat4x4_mult_vec(struct mat4x4 mat, struct vec3 vec)
{
	struct vec3 v;

	v.x = mat.vals[0] * vec.x + mat.vals[4] * vec.y + mat.vals[8] * vec.z + mat.vals[12];
	v.y = mat.vals[1] * vec.x + mat.vals[5] * vec.y + mat.vals[9] * vec.z + mat.vals[13];
	v.z = mat.vals[2] * vec.x + mat.vals[6] * vec.y + mat.vals[10] * vec.z + mat.vals[14];

	return v;
}

vec3 mat4x4_mult_vec3(struct mat4x4 mat, struct vec3 vec, double *last)
{
	struct vec3 v;

	v.x = mat.vals[0] * vec.x + mat.vals[4] * vec.y + mat.vals[8] * vec.z + mat.vals[12] * *last;
	v.y = mat.vals[1] * vec.x + mat.vals[5] * vec.y + mat.vals[9] * vec.z + mat.vals[13] * *last;
	v.z = mat.vals[2] * vec.x + mat.vals[6] * vec.y + mat.vals[10] * vec.z + mat.vals[14] * *last;
	
	*last = mat.vals[3] * vec.x + mat.vals[7] * vec.y
	+ mat.vals[11] * vec.z + mat.vals[15] * *last;
	
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

mat4x4 mat4x4_frustum(float left, float right, float top,
                      float bottom, float near, float far)
{
	mat4x4 mat = make_mat4x4();

	float r_width  = 1.0f / (right - left);
	float r_height = 1.0f / (top - bottom);
	float r_depth  = 1.0f / (far - near);
	float x =  2.0f * (r_width);
	float y =  2.0f * (r_height);
	float z =  2.0f * (r_depth);
	float A = (right + left) * r_width;
	float B = (top + bottom) * r_height;
	float C = (far + near) * r_depth;
//	mat.vals[0] = x;
//	mat.vals[3] = -A;
//	mat.vals[5] = y;
//	mat.vals[7] = -B;
//	mat.vals[10] = -z;
//	mat.vals[11] = -C;
	
	mat.vals[0] = 2 * near / (right - left);
	mat.vals[2] = (right + left) / (right - left);
	mat.vals[5] = 2 * near / (top - bottom);
	mat.vals[6] = (top + bottom) / (top - bottom);
	mat.vals[10] = - (far + near) / (far - near);
	mat.vals[11] = - 2 * far * near / (far - near);
	mat.vals[14] = -1;
	
	return mat;
}

mat4x4 mat4x4_ortho(float left, float right, float top,
                    float bottom, float near, float far)
{
	mat4x4 mat = make_mat4x4();

	mat.vals[0] = 2 / (right - left);
	mat.vals[3] = -(left + right) / (right - left);
	mat.vals[5] = 2 / (top - bottom);
	mat.vals[7] = -(top + bottom) / (top - bottom);
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

mat4x4 mat4x4_from_rot_trans(mat3x3 rot, vec3 trans)
{
	mat4x4 rot4 = make_mat4x4();

	rot4.vals[0] = rot.vals[0];
	rot4.vals[1] = rot.vals[1];
	rot4.vals[2] = rot.vals[2];

	rot4.vals[4] = rot.vals[3];
	rot4.vals[5] = rot.vals[4];
	rot4.vals[6] = rot.vals[5];

	rot4.vals[8] = rot.vals[6];
	rot4.vals[9] = rot.vals[7];
	rot4.vals[10] = rot.vals[8];
	
	mat4x4 trans4 = make_mat4x4();
	mat4x4_translate(&trans4, trans);

	return mat4x4_mult_mat4x4(trans4, rot4);
}

mat3x3 mat4x4_get_rot(mat4x4 &mat)
{
	mat3x3 rot = make_mat3x3();
	rot.vals[0] = mat.vals[0];
	rot.vals[1] = mat.vals[1];
	rot.vals[2] = mat.vals[2];
	rot.vals[3] = mat.vals[4];
	rot.vals[4] = mat.vals[5];
	rot.vals[5] = mat.vals[6];
	rot.vals[6] = mat.vals[8];
	rot.vals[7] = mat.vals[9];
	rot.vals[8] = mat.vals[10];
	return rot;
}

mat4x4 mat4x4_transpose(mat4x4 m)
{
	mat4x4 t;
	t.vals[0] = m.vals[0];
	t.vals[4] = m.vals[1];
	t.vals[8] = m.vals[2];
	t.vals[12] = m.vals[3];

	t.vals[1] = m.vals[4];
	t.vals[5] = m.vals[5];
	t.vals[9] = m.vals[6];
	t.vals[13] = m.vals[7];

	t.vals[2] = m.vals[8];
	t.vals[6] = m.vals[9];
	t.vals[10] = m.vals[10];
	t.vals[14] = m.vals[11];

	t.vals[3] = m.vals[12];
	t.vals[7] = m.vals[13];
	t.vals[11] = m.vals[14];
	t.vals[15] = m.vals[15];

	return t;
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

mat4x4 mat4x4_inverse(mat4x4 &mat)
{
	vec3 trans = make_vec3(mat.vals[12], mat.vals[13],
	                        mat.vals[14]);
	vec3_mult(&trans, -1);
	mat3x3 rot = mat4x4_get_rot(mat);	
	mat3x3 inv = mat3x3_transpose(rot);
	mat4x4 both = mat4x4_from_rot_trans(inv, trans);
	return both;
}
