// vagabond
// Copyright (C) 2019 Helen Ginn
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

#include "Blob.h"
#include "FFT.h"
#include "Options.h"
#include <cfloat>
#include <iostream>

Blob::Blob()
{
	_scale = 0.2;
	_random = rand() % 1000000;
	_added = false;
}

bool xPolygon(vec3 point, vec3 *vs)
{
	bool c = false;
	
	for (int i = 0, j = 2; i < 3; j = i++) 
	{
		if (((vs[i].z > point.z) != (vs[j].z > point.z))
		    && (point.y < (vs[j].y - vs[i].y) * (point.z - vs[i].z)
		    / (vs[j].z - vs[i].z) + vs[i].y)) 
		{
			c = !c;
		}
	}

	return c;
}

bool yPolygon(vec3 point, vec3 *vs)
{
	bool c = false;
	
	for (int i = 0, j = 2; i < 3; j = i++) 
	{
		if (((vs[i].z > point.z) != (vs[j].z > point.z))
		    && (point.x < (vs[j].x - vs[i].x) * (point.z - vs[i].z)
		    / (vs[j].z - vs[i].z) + vs[i].x)) 
		{
			c = !c;
		}
	}

	return c;
}

bool zPolygon(vec3 point, vec3 *vs)
{
	bool c = false;
	
	for (int i = 0, j = 2; i < 3; j = i++) 
	{
		if (((vs[i].y > point.y) != (vs[j].y > point.y))
		    && (point.x < (vs[j].x - vs[i].x) * (point.y - vs[i].y)
		    / (vs[j].y - vs[i].y) + vs[i].x)) 
		{
			c = !c;
		}
	}

	return c;
}

bool Blob::polygonIncludes(vec3 point, unsigned long *trio)
{
	vec3 vs[3];
	vs[0] = _vertices[trio[0]]; 
	vs[1] = _vertices[trio[1]];
	vs[2] = _vertices[trio[2]];

	double xmin = std::min(std::min(vs[0].x, vs[1].x), vs[2].x);
	double ymin = std::min(std::min(vs[0].y, vs[1].y), vs[2].y);
	double zmin = std::min(std::min(vs[0].z, vs[1].z), vs[2].z);
	
	double xmax = std::max(std::max(vs[0].x, vs[1].x), vs[2].x);
	double ymax = std::max(std::max(vs[0].y, vs[1].y), vs[2].y);
	double zmax = std::max(std::max(vs[0].z, vs[1].z), vs[2].z);
	
	double zdiff = zmax - zmin;
	double ydiff = ymax - ymin;
	double xdiff = xmax - xmin;
	
	if (zdiff < ydiff && zdiff < xdiff)
	{
		return zPolygon(point, (vec3 *)vs);
	}
	else if (ydiff < zdiff && ydiff < xdiff)
	{
		return yPolygon(point, (vec3 *)vs);
	}
	else
	{
		return xPolygon(point, (vec3 *)vs);
	}
}

vec3 Blob::rayTraceToPlane(vec3 point, unsigned long *trio, vec3 dir, bool *b)
{
	vec3 vs[3];
	vs[0] = _vertices[trio[0]]; 
	vs[1] = _vertices[trio[1]];
	vs[2] = _vertices[trio[2]];

	vec3 diff1 = vec3_subtract_vec3(vs[1], vs[0]);
	vec3 diff2 = vec3_subtract_vec3(vs[2], vs[0]);
	vec3 cross = vec3_cross_vec3(diff1, diff2);
	vec3_set_length(&cross, 1); 
	
	double denom = vec3_dot_vec3(dir, cross);
	vec3 subtract = vec3_subtract_vec3(vs[0], point);
	double nom = vec3_dot_vec3(subtract, cross);
	double d = nom / denom;
	
	vec3_mult(&dir, d);
	vec3_add_to_vec3(&point, dir);
	
	*b = (d < 0);
	return point;
}

bool Blob::pointInside(vec3 point)
{
	vec3 dir = make_vec3(0, 0, 1);
	bool c = false;

	for (size_t i = 0; i < _indices.size(); i += 3)
	{
		unsigned long *ptr = &_indices[i];
		bool backwards = false;
		
		if ((point.x < _vertices[_indices[i]].x &&
		    point.x < _vertices[_indices[i+1]].x &&
		     point.x < _vertices[_indices[i+2]].x) ||
		(point.x > _vertices[_indices[i]].x &&
		 point.x > _vertices[_indices[i+1]].x &&
		 point.x > _vertices[_indices[i+2]].x))
		{
			continue;
		}

		if ((point.y < _vertices[_indices[i]].y &&
		    point.y < _vertices[_indices[i+1]].y &&
		     point.y < _vertices[_indices[i+2]].y) ||
		(point.y > _vertices[_indices[i]].y &&
		 point.y > _vertices[_indices[i+1]].y &&
		 point.y > _vertices[_indices[i+2]].y))
		{
			continue;
		}

		vec3 ray = rayTraceToPlane(point, ptr, dir, &backwards);

		if (backwards)
		{
			continue;
		}

		vec3 vs[3];
		vs[0] = _vertices[ptr[0]];
		vs[1] = _vertices[ptr[1]];
		vs[2] = _vertices[ptr[2]];
		bool inside = zPolygon(ray, (vec3 *)vs);

		if (inside)
		{
			c = !c;
		}
	}

	return c;
}


vec3 Blob::closestRayTraceToPlane(vec3 point, unsigned long *trio)
{
	vec3 vs[3];
	vs[0] = _vertices[trio[0]]; 
	vs[1] = _vertices[trio[1]];
	vs[2] = _vertices[trio[2]];

	vec3 diff1 = vec3_subtract_vec3(vs[1], vs[0]);
	vec3 diff2 = vec3_subtract_vec3(vs[2], vs[0]);
	vec3 cross = vec3_cross_vec3(diff1, diff2);
	vec3_set_length(&cross, 1); 
	
	vec3 subtract = vec3_subtract_vec3(vs[0], point);
	double d = vec3_dot_vec3(subtract, cross);
	
	vec3_mult(&cross, d);
	vec3_add_to_vec3(&point, cross);
	
	return point;
}

vec3 Blob::findClosestVecToSurface(vec3 &interior)
{
	double distance = FLT_MAX;
	vec3 closest = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);

	for (size_t i = 0; i < _indices.size(); i += 3)
	{
		vec3 close = closestRayTraceToPlane(interior, &_indices[i]);
		
		bool included = polygonIncludes(close, &_indices[i]);
		
		if (!included)
		{
			continue;
		}

		vec3_subtract_from_vec3(&close, interior);

		double l = vec3_length(close);
		
		if (l >= distance || l != l)
		{
			continue;
		}
		
		closest = close;
		distance = l;
	}

	vec3_mult(&closest, -1);
	return closest;
}

void Blob::getLoopBoundaries(vec3 *min, vec3 *max)
{
	for (size_t i = 0; i < _vertices.size(); i++)
	{
		if (_vertices[i].x < min->x) { min->x = _vertices[i].x; }
		if (_vertices[i].x > max->x) { max->x = _vertices[i].x; }

		if (_vertices[i].y < min->y) { min->y = _vertices[i].y; }
		if (_vertices[i].y > max->y) { max->y = _vertices[i].y; }

		if (_vertices[i].z < min->z) { min->z = _vertices[i].z; }
		if (_vertices[i].z > max->z) { max->z = _vertices[i].z; }
	}
}

void Blob::estimateScale(VagFFTPtr fft)
{
	vec3 min, max;
	min = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	max = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	getLoopBoundaries(&min, &max);

	min.x -= 2; min.y -= 2; min.z -= 2;
	max.x += 2; max.y += 2; max.z += 2;

	double sum = 0;
	double count = 0;

	mat3x3 recip = fft->getRecipBasis();
	double cubeDim = Options::getActiveCrystal()->getProteinSampling();
	
	for (double z = min.z; z < max.z; z += cubeDim)
	{
		for (double y = min.y; y < max.y; y += cubeDim)
		{
			for (double x = min.x; x < max.x; x += cubeDim)
			{
				vec3 point = make_vec3(x, y, z);

				if (!pointInside(point))
				{
					continue;
				}
				
				mat3x3_mult_vec(recip, &point);
				double density = fft->getCompFromFrac(point, 0);

				sum += density;
				count += 1;
			}
		}
	}
	
	_scale = sum / count;
	_scale /= 4;
	_scale = 0.2;
	std::cout << "Scale for blob: " << _scale << std::endl;
}

void Blob::addToCubicMap(VagFFTPtr fft)
{
	vec3 min, max;
	min = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	max = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	getLoopBoundaries(&min, &max);

	min.x -= 2; min.y -= 2; min.z -= 2;
	max.x += 2; max.y += 2; max.z += 2;
	
	double plateau = 2;
	
	mat3x3 recip = fft->getRecipBasis();
	double cubeDim = Options::getActiveCrystal()->getProteinSampling();
	vec3 origin = fft->origin();
	
	for (double z = min.z; z < max.z; z += cubeDim)
	{
		for (double y = min.y; y < max.y; y += cubeDim)
		{
			for (double x = min.x; x < max.x; x += cubeDim)
			{
				vec3 point = make_vec3(x, y, z);

				double density = 1;
				
				if (!pointInside(point))
				{
					continue;
				}
				else
				{
					vec3 toSurface = findClosestVecToSurface(point);

					double dist = vec3_length(toSurface);
					if (dist < plateau)
					{
						density = dist / plateau;
					}
				}
				
				density *= _scale;
				vec3_subtract_from_vec3(&point, origin);

				mat3x3_mult_vec(recip, &point);
				fft->addInterpolatedToReal(point.x, point.y, point.z, 
				                           density, -1);
			}
		}
	}
}

std::string Blob::getParserIdentifier()
{
	return "Blob_" + i_to_str(_random);
}

void Blob::addProperties()
{
	addVec3ArrayProperty("vertices", &_vertices);
	addIntArrayProperty("indices", &_indices);
	addDoubleProperty("scale", &_scale);
	addIntProperty("random", &_random);

}

void Blob::postParseTidy()
{
	if (_added)
	{
		return;
	}

	Options::loadBlobInGUI(this);
	_added = true;
}
