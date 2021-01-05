//
//  mat3x3.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include <hcsrc/vec3.h>
#include <hcsrc/mat3x3.h>
#include "../libccp4/csymlib.h"

/** mat3x3 matrix from CCP4 symmetry operation as a 
 * rotation matrix */
inline mat3x3 mat3x3_from_ccp4(CSym::ccp4_symop symop)
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


