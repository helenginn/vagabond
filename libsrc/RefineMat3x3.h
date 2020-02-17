// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#ifndef __vagabond__RefineMat3x3__
#define __vagabond__RefineMat3x3__

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
	
	void setMat3x3(mat3x3 mat)
	{
		_mat = mat;
	}
	
	static RefineMat3x3 *toMat(void *object)
	{
		return static_cast<RefineMat3x3 *>(object);
	}

	static void setTensor11(void *object, double value)
	{
		toMat(object)->_mat.vals[0] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor12(void *object, double value)
	{
		toMat(object)->_mat.vals[1] = value;
		toMat(object)->_mat.vals[3] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor21(void *object, double value)
	{
		toMat(object)->_mat.vals[3] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor13(void *object, double value)
	{
		toMat(object)->_mat.vals[2] = value;
		toMat(object)->_mat.vals[6] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor31(void *object, double value)
	{
		toMat(object)->_mat.vals[6] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor22(void *object, double value)
	{
		toMat(object)->_mat.vals[4] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor32(void *object, double value)
	{
		toMat(object)->_mat.vals[7] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor23(void *object, double value)
	{
		toMat(object)->_mat.vals[5] = value;
		toMat(object)->_mat.vals[7] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setTensor33(void *object, double value)
	{
		toMat(object)->_mat.vals[8] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue11(void *object, double value)
	{
		toMat(object)->_mat.vals[0] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue12(void *object, double value)
	{
		toMat(object)->_mat.vals[1] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue13(void *object, double value)
	{
		toMat(object)->_mat.vals[2] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue21(void *object, double value)
	{
		toMat(object)->_mat.vals[3] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue22(void *object, double value)
	{
		toMat(object)->_mat.vals[4] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue23(void *object, double value)
	{
		toMat(object)->_mat.vals[5] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue31(void *object, double value)
	{
		toMat(object)->_mat.vals[6] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue32(void *object, double value)
	{
		toMat(object)->_mat.vals[7] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
	}

	static void setValue33(void *object, double value)
	{
		toMat(object)->_mat.vals[8] = value;
		if (*toMat(object)->_clean)
		{
			(*toMat(object)->_clean)(toMat(object)->_parent);
		}
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
	
	static double getValue11(void *object)
	{
		return toMat(object)->_mat.vals[0];
	}

	static double getValue12(void *object)
	{
		return toMat(object)->_mat.vals[1];
	}

	static double getValue13(void *object)
	{
		return toMat(object)->_mat.vals[2];
	}

	static double getValue21(void *object)
	{
		return toMat(object)->_mat.vals[3];
	}

	static double getValue22(void *object)
	{
		return toMat(object)->_mat.vals[4];
	}

	static double getValue23(void *object)
	{
		return toMat(object)->_mat.vals[5];
	}

	static double getValue31(void *object)
	{
		return toMat(object)->_mat.vals[6];
	}

	static double getValue32(void *object)
	{
		return toMat(object)->_mat.vals[7];
	}

	static double getValue33(void *object)
	{
		return toMat(object)->_mat.vals[8];
	}
	
	void setZero()
	{
		mat3x3_mult_scalar(&_mat, 0);
	}
	
	void addMajorAxesToStrategy(RefinementStrategyPtr strategy, double step,
	                            double tol, std::string prefix)
	{
		strategy->addParameter(this, getTensor11, setTensor11, step,
		                       tol, prefix + "_t11");
		strategy->addParameter(this, getTensor22, setTensor22, step,
		                       tol, prefix + "_t22");
		strategy->addParameter(this, getTensor33, setTensor33, step,
		                       tol, prefix + "_t33");

	}

	void addTensorToStrategy(RefinementStrategyPtr strategy, double step,
	                         double tol, std::string prefix)
	{
		strategy->addParameter(this, getTensor11, setTensor11, step,
		                       tol, prefix + "_t11");
		strategy->addParameter(this, getTensor22, setTensor22, step,
		                       tol, prefix + "_t22");
		strategy->addParameter(this, getTensor33, setTensor33, step,
		                       tol, prefix + "_t33");
		strategy->addParameter(this, getTensor12, setTensor12, step,
		                       tol, prefix + "_t12");
		strategy->addParameter(this, getTensor13, setTensor13, step,
		                       tol, prefix + "_t13");
		strategy->addParameter(this, getTensor23, setTensor23, step,
		                       tol, prefix + "_t23");
	}
	
	void addMatrixToStrategy(RefinementStrategyPtr strategy, double step,
	                   double tol, std::string prefix)
	{
		strategy->addParameter(this, getValue11, setValue11, step,
		                       tol, prefix + "_t11");
		strategy->addParameter(this, getValue12, setValue12, step,
		                       tol, prefix + "_t12");
		strategy->addParameter(this, getValue13, setValue13, step,
		                       tol, prefix + "_t13");
		strategy->addParameter(this, getValue21, setValue21, step,
		                       tol, prefix + "_t21");
		strategy->addParameter(this, getValue22, setValue22, step,
		                       tol, prefix + "_t22");
		strategy->addParameter(this, getValue23, setValue23, step,
		                       tol, prefix + "_t23");
		strategy->addParameter(this, getValue31, setValue31, step,
		                       tol, prefix + "_t31");
		strategy->addParameter(this, getValue32, setValue32, step,
		                       tol, prefix + "_t32");
		strategy->addParameter(this, getValue33, setValue33, step,
		                       tol, prefix + "_t33");

	}
	
private:
	mat3x3 _mat;
	void *_parent;
	CleanUp _clean;
};

#endif

