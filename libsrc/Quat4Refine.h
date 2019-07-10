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

#ifndef __Vagabond__Quat4Refine__
#define __Vagabond__Quat4Refine__

#include "quat4.h"
#include "RefinementStrategy.h"

class Quat4Refine
{
public:
	Quat4Refine(void *parent = NULL)
	{
		_parent = parent;
		_quat = empty_quat4();
	}
	
	quat4 getQuat4()
	{
		return _quat;
	}
	
	vec3 getVec3()
	{
		return make_vec3(_quat.x, _quat.y, _quat.z);
	}
	
	static Quat4Refine *toQuat(void *object)
	{
		return static_cast<Quat4Refine *>(object);
	}
	
	static void setX(void *object, double value)
	{
		toQuat(object)->_quat.x = value;
	}
	
	static void setY(void *object, double value)
	{
		toQuat(object)->_quat.y = value;
	}
	
	static void setZ(void *object, double value)
	{
		toQuat(object)->_quat.z = value;
	}
	
	static void setT(void *object, double value)
	{
		toQuat(object)->_quat.t = value;
	}
	
	static double getX(void *object)
	{
		return toQuat(object)->_quat.x;
	}
	
	static double getY(void *object)
	{
		return toQuat(object)->_quat.y;
	}
	
	static double getZ(void *object)
	{
		return toQuat(object)->_quat.z;
	}
	
	static double getT(void *object)
	{
		return toQuat(object)->_quat.t;
	}

	void addQuatToStrategy(RefinementStrategyPtr strategy, double step,
	                       double tol, std::string prefix)
	{
		strategy->addParameter(this, getX, setX, step, tol, prefix + "_x");
		strategy->addParameter(this, getY, setY, step, tol, prefix + "_y");
		strategy->addParameter(this, getZ, setZ, step, tol, prefix + "_z");
		strategy->addParameter(this, getT, setT, step, tol, prefix + "_t");
	}

	void addVecToStrategy(RefinementStrategyPtr strategy, double step,
	                       double tol, std::string prefix)
	{
		strategy->addParameter(this, getX, setX, step, tol, prefix + "_x");
		strategy->addParameter(this, getY, setY, step, tol, prefix + "_y");
		strategy->addParameter(this, getZ, setZ, step, tol, prefix + "_z");
	}

	void addVec2ToStrategy(RefinementStrategyPtr strategy, double step,
	                       double tol, std::string prefix)
	{
		strategy->addParameter(this, getX, setX, step, tol, prefix + "_x");
		strategy->addParameter(this, getY, setY, step, tol, prefix + "_y");
	}
private:
	void *_parent;
	quat4 _quat;
};

#endif
