//
//  RefineMat3x3.h
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "mat3x3.h"
#include "RefinementStrategy.h"

typedef void (*CleanUp)(void *);


class RefineMat3x3
{
public:
	RefineMat3x3(void *parent, CleanUp clean)
	{
		_parent = parent;
		_clean = clean;
		_mat = make_mat3x3();
	}
	
	mat3x3 getMat3x3()
	{
		return _mat;
	}
	
	mat3x3 *getMat3x3Ptr()
	{
		return &_mat;
	}
	
	static RefineMat3x3 *toMat(void *object)
	{
		return static_cast<RefineMat3x3 *>(object);
	}

	static void setTensor11(void *object, double value)
	{
		toMat(object)->_mat.vals[0] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor12(void *object, double value)
	{
		toMat(object)->_mat.vals[1] = value;
		toMat(object)->_mat.vals[3] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor21(void *object, double value)
	{
		toMat(object)->_mat.vals[3] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor13(void *object, double value)
	{
		toMat(object)->_mat.vals[2] = value;
		toMat(object)->_mat.vals[6] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor31(void *object, double value)
	{
		toMat(object)->_mat.vals[6] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor22(void *object, double value)
	{
		toMat(object)->_mat.vals[4] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor32(void *object, double value)
	{
		toMat(object)->_mat.vals[7] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor23(void *object, double value)
	{
		toMat(object)->_mat.vals[5] = value;
		toMat(object)->_mat.vals[7] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static void setTensor33(void *object, double value)
	{
		toMat(object)->_mat.vals[8] = value;
		(*toMat(object)->_clean)(toMat(object)->_parent);
	}

	static double getTensor11(void *object)
	{
		return toMat(object)->_mat.vals[0];
	}

	static double getTensor21(void *object)
	{
		return toMat(object)->_mat.vals[3];
	}

	static double getTensor12(void *object)
	{
		return toMat(object)->_mat.vals[1];
	}

	static double getTensor31(void *object)
	{
		return toMat(object)->_mat.vals[6];
	}

	static double getTensor13(void *object)
	{
		return toMat(object)->_mat.vals[2];
	}

	static double getTensor22(void *object)
	{
		return toMat(object)->_mat.vals[4];
	}

	static double getTensor23(void *object)
	{
		return toMat(object)->_mat.vals[5];
	}

	static double getTensor32(void *object)
	{
		return toMat(object)->_mat.vals[7];
	}

	static double getTensor33(void *object)
	{
		return toMat(object)->_mat.vals[8];
	}
	
	void setZero()
	{
		mat3x3_mult_scalar(&_mat, 0);
	}

	void addToStrategy(RefinementStrategyPtr strategy, double step,
	                   double tol, std::string prefix)
	{
		strategy->addParameter(this, getTensor11, setTensor11, step,
		                       tol, prefix + "_t11");
		strategy->addParameter(this, getTensor12, setTensor12, step,
		                       tol, prefix + "_t12");
		strategy->addParameter(this, getTensor13, setTensor13, step,
		                       tol, prefix + "_t13");
		strategy->addParameter(this, getTensor22, setTensor22, step,
		                       tol, prefix + "_t22");
		strategy->addParameter(this, getTensor23, setTensor23, step,
		                       tol, prefix + "_t23");
		strategy->addParameter(this, getTensor33, setTensor33, step,
		                       tol, prefix + "_t33");
	}
private:
	mat3x3 _mat;
	void *_parent;
	CleanUp _clean;
};
