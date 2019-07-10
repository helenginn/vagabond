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

#ifndef __Vagabond__quat4__
#define __Vagabond__quat4__

#include "vec3.h"

struct quat4
{
	double t;
	double x;
	double y;
	double z;
};

inline quat4 empty_quat4()
{
	struct quat4 quat;
	quat.x = 0;
	quat.y = 0;
	quat.z = 0;
	quat.t = 1;
	
	return quat;
}

inline quat4 make_quat4(double x, double y, double z, double t)
{
	struct quat4 quat;
	quat.x = x;
	quat.y = y;
	quat.z = z;
	quat.t = t;

	return quat;
}

inline double quat4_angle(quat4 &quat)
{
	double t = quat.t;
	double ang = acos(t) * 2;
	
	return ang;
}

inline vec3 quat4_axis(quat4 &quat)
{
	vec3 vec = make_vec3(quat.x, quat.y, quat.z);
	vec3_set_length(&vec, 1);
	
	return vec;
}

inline quat4 quat4_mult_quat4(quat4 &q1, quat4 &q2)
{
	struct quat4 prod;
	prod.t = q1.t * q2.t - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
	prod.x = q1.t * q2.x + q1.x * q2.t + q1.y * q2.z - q1.z * q2.y;
	prod.y = q1.t * q2.y + q1.x * q2.z + q1.y * q2.t - q1.z * q2.x;
	prod.y = q1.t * q2.z + q1.x * q2.y + q1.y * q2.x - q1.z * q2.t;
	
	return prod;
}

inline quat4 quat4_inverse(quat4 &q)
{
	quat4 inv = make_quat4(-q.x, -q.y, -q.z, q.t);
	return inv;
}

inline vec3 quat4_rot_vec3(quat4 &q, vec3 &v)
{
	quat4 qi = quat4_inverse(q);
	quat4 qv = make_quat4(v.x, v.y, v.z, 0);
	
	quat4 first = quat4_mult_quat4(qv, qi);
	quat4 prod  = quat4_mult_quat4(q, first);
	
	return make_vec3(prod.x, prod.y, prod.z);
}

#endif
