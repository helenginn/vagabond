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
#include "RefinementStepSearch.h"
#include "RefinementNelderMead.h"
#include <iostream>

void Kabsch::setAtoms(std::vector<vec3> aAtoms, std::vector<vec3> bAtoms)
{
	if (aAtoms.size() != bAtoms.size())
	{
		shout_at_helen("Attempting Kabsch algorithm on differing atom counts!");
	}

	_vecs[0] = aAtoms;
	_vecs[1] = bAtoms;

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

void Kabsch::setWeights(std::vector<double> &weights)
{
	_weights = weights;
}

vec3 Kabsch::fixCentroids()
{
	vec3 aSum = make_vec3(0, 0, 0);
	vec3 bSum = make_vec3(0, 0, 0);

	for (int i = 0; i < _positions[0].size(); i++)
	{
		aSum = vec3_add_vec3(aSum, _vecs[0][i]);
		bSum = vec3_add_vec3(bSum, _vecs[1][i]);
	}

	double mult = 1 / (double)_positions[0].size();
	vec3_mult(&aSum, mult);
	vec3_mult(&bSum, mult);

	for (int i = 0; i < _positions[0].size(); i++)
	{
		_positions[0][i][0] -= aSum.x;
		_positions[0][i][1] -= aSum.y;
		_positions[0][i][2] -= aSum.z;
		_positions[1][i][0] -= bSum.x;
		_positions[1][i][1] -= bSum.y;
		_positions[1][i][2] -= bSum.z;

		_vecs[0][i] = vec3_subtract_vec3(_vecs[0][i], aSum);
		_vecs[1][i] = vec3_subtract_vec3(_vecs[1][i], bSum);
	}

	return aSum;
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

	double weightSum = 0;

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int k = 0; k < _positions[0].size(); k++)
			{
				double add = _positions[0][k][i] * _positions[1][k][j];
				double weight = 1;
				if (_weights.size() > k)
				{
					weight = _weights[i] * _weights[j];
				}
				add *= weight;
				weightSum += weight;
				matrix[i][j] += add;
			}
		}
	}

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			matrix[i][j] /= weightSum;
		}
	}

	_covariance = mat3x3_from_2d_array(matrix);

	vect w = (double *)malloc(sizeof(double) * 3);
	int success = svdcmp(matrix, 3, 3, w, v);

	if (!success)
	{
		_failed = true;
		_rotation = make_mat3x3();
		return _rotation;
	}

	mat3x3 myU = mat3x3_from_2d_array(matrix); // has been changed by svdcmp.
	mat3x3 myV = mat3x3_from_2d_array(v);

	mat3x3 vT = mat3x3_transpose(myV);
	mat3x3 mult = mat3x3_mult_mat3x3(myU, vT);
	double det = mat3x3_determinant(mult);

	if (det > 0)
	{
		_rotation = mult;

		free_2d_array(matrix);
		free_2d_array(v);
		free(w);

		return mult;
	}

	mat3x3 detMat = make_mat3x3();
	detMat.vals[8] = det;

	mat3x3 mat1 = mat3x3_mult_mat3x3(detMat, vT);
	/*** BUG HERE ***/
	mat3x3 mat2 = mat3x3_mult_mat3x3(myU, mat1);

	_rotation = mat2;

	free_2d_array(matrix);
	free_2d_array(v);
	free(w);
	
	return mat2;
}

mat3x3 Kabsch::findFinalTransform()
{
	for (int i = 0; i < _vecs[0].size(); i++)
	{
		mat3x3_mult_vec(_rotation, &_vecs[0][i]);
	}

	const double dims[] = {1, 1, 1, 90, 90, 90};
	_unitCell = std::vector<double>(dims, dims+6);

	// now as close as possible without doing a funny scaling thing.

	RefinementStepSearchPtr mead = RefinementStepSearchPtr(new RefinementStepSearch());
	//NelderMeadPtr mead = NelderMeadPtr(new NelderMead());
	mead->setCycles(50);
	mead->addParameter(this, getUnitCellA, setUnitCellA, 0.02, 0.0001);
	mead->addParameter(this, getUnitCellB, setUnitCellB, 0.02, 0.0001);
	mead->addParameter(this, getUnitCellC, setUnitCellC, 0.02, 0.0001);
	mead->addParameter(this, getUnitCellAlpha, setUnitCellAlpha, 0.2, 0.0001);
	mead->addParameter(this, getUnitCellBeta, setUnitCellBeta, 0.2, 0.0001);
	mead->addParameter(this, getUnitCellGamma, setUnitCellGamma, 0.2, 0.0001);
	mead->setEvaluationFunction(Kabsch::score, this);

	mead->refine();

	return mat3x3_from_unit_cell(&(_unitCell[0]));
}

double Kabsch::score(void *object)
{
	Kabsch *me = static_cast<Kabsch *>(object);

	mat3x3 unit_cell = mat3x3_from_unit_cell(&(me->_unitCell[0]));

	double diffSum = 0;

	for (int i = 0; i < me->_vecs[0].size(); i++)
	{
		vec3 tmp = me->_vecs[0][i];
		mat3x3_mult_vec(unit_cell, &tmp);

		vec3 difference = vec3_subtract_vec3(me->_vecs[1][i], tmp);
		diffSum += vec3_sqlength(difference);
	}

	diffSum /= me->_vecs[0].size();

	return diffSum;
}
