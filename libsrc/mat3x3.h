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

struct mat3x3
{
	double vals[9];
};


mat3x3 mat3x3_inverse(mat3x3 &mat);
mat3x3 mat3x3_from_unit_cell(double a, double b, double c, double alpha, double beta, double gamma);
mat3x3 make_mat3x3();
void mat3x3_mult_vec(struct mat3x3 mat, struct vec3 *vec);
void mat3x3_scale(mat3x3 *mat, double a, double b, double c);
double mat3x3_length(mat3x3 &mat, int index);
mat3x3 mat3x3_transpose(mat3x3 &mat);
double mat3x3_determinant(mat3x3 &mat);

#endif /* defined(__vagabond__mat3x3__) */
