//
//  mat4x4.hpp
//  VagabondViewer
//
//  Created by Helen Ginn on 02/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef mat4x4_hpp
#define mat4x4_hpp

#include <stdio.h>
#include <string>
#include "vec3.h"
#include "mat3x3.h"

struct mat4x4
{
	float vals[16];
};

std::string mat4x4_desc(mat4x4 mat);
mat4x4 make_mat4x4();
mat4x4 mat4x4_frustum(float left, float right, float top,
                      float bottom, float near, float far);
mat4x4 mat4x4_ortho(float left, float right, float top,
                    float bottom, float near, float far);

void mat4x4_rotate(mat4x4 *mat, double alpha, double beta, double gamma);
void mat4x4_translate(mat4x4 *mat, vec3 centre);
mat4x4 mat4x4_from_rot_trans(mat3x3 rot, vec3 trans);
mat3x3 mat4x4_get_rot(mat4x4 &mat);
mat4x4 mat4x4_mult_mat4x4(mat4x4 l, mat4x4 r);
vec3 mat4x4_mult_vec(struct mat4x4 mat, struct vec3 vec);
void mat4x4_mult_scalar(mat4x4 *mat, double scalar);
vec3 mat4x4_mult_vec3(struct mat4x4 mat, struct vec3 vec, 
                      double *last);
mat4x4 mat4x4_inverse(mat4x4 &mat);
mat4x4 mat4x4_transpose(mat4x4 m);

#endif /* mat4x4_hpp */
