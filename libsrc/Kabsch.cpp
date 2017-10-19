//
//  Kabsch.cpp
//  vagabond
//
//  Created by Helen Ginn on 18/10/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "Kabsch.h"
#include "Shouter.h"
#include "Model.h"
#include "matrix.h"
#include "svdcmp.h"

void Kabsch::setAtoms(std::vector<vec3> aAtoms, std::vector<vec3> bAtoms)
{
	if (aAtoms.size() != bAtoms.size())
	{
		shout_at_helen("Attempting Kabsch algorithm on differing atom counts!");
	}

	for (int i = 0; i < aAtoms.size(); i++)
	{
		vec3 aPosVec = aAtoms[i];
		vec3 bPosVec = bAtoms[i];

		std::vector<double> aPos, bPos;
		aPos.push_back(aPosVec.x);
		aPos.push_back(aPosVec.y);
		aPos.push_back(aPosVec.z);
		bPos.push_back(bPosVec.x);
		bPos.push_back(bPosVec.y);
		bPos.push_back(bPosVec.z);

		_positions[0].push_back(aPos);
		_positions[1].push_back(bPos);
	}
}

mat3x3 Kabsch::run()
{
	if (!_positions[0].size())
	{
		shout_at_helen("Kabsch run with no atoms");
	}

	mat matrix = (double **)malloc(sizeof(double *) * 3);
	mat v = (double **)malloc(sizeof(double *) * 3);

	for (int i = 0; i < 3; i++)
	{
		matrix[i] = (double *)malloc(sizeof(double) * 3);
		v[i] = (double *)malloc(sizeof(double) * 3);
		memset(matrix[i], 0, sizeof(double) * 3);
		memset(v[i], 0, sizeof(double) * 3);
	}

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int k = 0; k < _positions[0].size(); k++)
			{
				matrix[i][j] += _positions[0][k][i] * _positions[1][k][j];
			}
		}
	}

	vect w = (double *)malloc(sizeof(double) * 3);
	svdcmp(matrix, 3, 3, w, v);

	mat3x3 myU = make_mat3x3();
	mat3x3 myV = make_mat3x3();

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			myU.vals[j * 3 + i] = matrix[i][j];
			myV.vals[j * 3 + i] = v[i][j];
		}
	}

	mat3x3 vT = mat3x3_transpose(myV);
	mat3x3 mult = mat3x3_mult_mat3x3(myU, vT);
	double det = mat3x3_determinant(mult);

	mat3x3 detMat = make_mat3x3();
	detMat.vals[8] = det;

	mat3x3 mat1 = mat3x3_mult_mat3x3(detMat, vT);
	mat3x3 mat2 = mat3x3_mult_mat3x3(myU, mat1);

	return mat2;
}
