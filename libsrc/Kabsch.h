//
//  Kabsch.hpp
//  vagabond
//
//  Created by Helen Ginn on 18/10/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef Kabsch_hpp
#define Kabsch_hpp

#include <stdio.h>
#include <vector>
#include "vec3.h"
#include "shared_ptrs.h"
#include "mat3x3.h"

/**
 * \class Kabsch
 * \brief Takes arrays of a series of vec3 points and generates
 * rotations/translations in order to superimpose the arrays.
 *
 * Used primarily in order to superimpose polymers on top of each other.
 * Named after Wolfgang Kabsch's algorithm to superimpose proteins in a
 * single step using SVD. This implementation superimposes them against the
 * average structure, which means that theoretically this needs to be
 * repeated several times to converge. In practice, the second cycle turns
 * out to be negligible.
 */


class Kabsch
{
public:
	Kabsch()
	{
		_covariance = make_mat3x3();
		_rotation = make_mat3x3();
		_failed = false;
	}

	void setWeights(std::vector<double> &weights);
	void setAtoms(std::vector<vec3> aAtoms, std::vector<vec3> bAtoms);
	mat3x3 run();
	mat3x3 findFinalTransform();
	vec3 fixCentroids();

	mat3x3 getCovariance()
	{
		return _covariance;
	}

	bool didFail()
	{
		return _failed;
	}

	static double score(void *object);

	static void setUnitCellA(void *object, double a)
	{
		static_cast<Kabsch *>(object)->_unitCell[0] = a;
	}

	static void setUnitCellB(void *object, double b)
	{
		static_cast<Kabsch *>(object)->_unitCell[1] = b;
	}

	static void setUnitCellC(void *object, double c)
	{
		static_cast<Kabsch *>(object)->_unitCell[2] = c;
	}

	static void setUnitCellAlpha(void *object, double alpha)
	{
		static_cast<Kabsch *>(object)->_unitCell[3] = alpha;
	}

	static void setUnitCellBeta(void *object, double beta)
	{
		static_cast<Kabsch *>(object)->_unitCell[4] = beta;
	}

	static void setUnitCellGamma(void *object, double gamma)
	{
		static_cast<Kabsch *>(object)->_unitCell[5] = gamma;
	}

	static double getUnitCellA(void *object)
	{
		return static_cast<Kabsch *>(object)->_unitCell[0];
	}

	static double getUnitCellB(void *object)
	{
		return static_cast<Kabsch *>(object)->_unitCell[1];
	}

	static double getUnitCellC(void *object)
	{
		return static_cast<Kabsch *>(object)->_unitCell[2];
	}

	static double getUnitCellAlpha(void *object)
	{
		return static_cast<Kabsch *>(object)->_unitCell[3];
	}

	static double getUnitCellBeta(void *object)
	{
		return static_cast<Kabsch *>(object)->_unitCell[4];
	}

	static double getUnitCellGamma(void *object)
	{
		return static_cast<Kabsch *>(object)->_unitCell[5];
	}

private:
	std::vector<std::vector<double> > _positions[2];
	std::vector<vec3> _vecs[2];
	std::vector<double> _weights;
	bool _failed;

	std::vector<double> _unitCell;
	mat3x3 _covariance;
	mat3x3 _rotation;
};

#endif /* Kabsch_hpp */
