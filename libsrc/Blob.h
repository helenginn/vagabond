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

#ifndef __vagabond__Blob__
#define __vagabond__Blob__

#include <hcsrc/vec3.h>
#include "Parser.h"
#include "shared_ptrs.h"

class Blob : public Parser
{
public:
	BlobPtr shared_from_this()
	{
		return ToBlobPtr(Parser::shared_from_this());
	}

	Blob();
	
	void addVertex(vec3 v)
	{
		_vertices.push_back(v);
	}
	
	size_t indexCount()
	{
		return _indices.size();
	}
	
	size_t vertexCount()
	{
		return _vertices.size();
	}
	
	unsigned long index(int i)
	{
		return _indices[i];
	}
	
	vec3 vertex(int i)
	{
		return _vertices[i];
	}
	
	void setAdded(bool added)
	{
		_added = true;
	}

	void addIndex(unsigned long idx)
	{
		_indices.push_back(idx);
	}

	void clear()
	{
		_vertices.clear();
		_indices.clear();
	}
	
	void multiplyScale(double mult)
	{
		_scale *= mult;
	}

	void addToCubicMap(VagFFTPtr fft);
	void getLoopBoundaries(vec3 *min, vec3 *max);
	void estimateScale(VagFFTPtr fft);
	
	std::string getClassName()
	{
		return "Blob";
	}
	
protected:
	virtual void addProperties();
	virtual void postParseTidy();
	std::string getParserIdentifier();
private:
	vec3 findClosestVecToSurface(vec3 &interior);
	vec3 rayTraceToPlane(vec3 point, unsigned long *trio, vec3 dir, bool *b);
	vec3 closestRayTraceToPlane(vec3 point, unsigned long *trio);
	bool polygonIncludes(vec3 point, unsigned long *trio);
	bool pointInside(vec3 point);
	bool _added;

	std::vector<vec3> _vertices;
	std::vector<unsigned long> _indices;

	double _scale;
	int _random;
};

#endif
