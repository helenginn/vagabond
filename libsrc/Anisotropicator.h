//
//  Anisotropicator.hpp
//  vagabond
//
//  Created by Helen Ginn on 28/10/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef Anisotropicator_hpp
#define Anisotropicator_hpp

#include <hcsrc/mat3x3.h>
#include <vector>
#include <stdio.h>

/**
 * \class Anisotropicator
 * \brief Class takes a bunch of points and performs SVD to recalculate a
 * tensor from them.
 */

class Anisotropicator
{
public:
	Anisotropicator();

	/* Should update the principle axes */
	void setTensor(mat3x3 tensor);
	/* Should update the tensor too */
	void setPoints(std::vector<vec3> points);

	vec3 getAxis(int i)
	{
		return _axes[i];
	}

	/** Return matrix describing the principle axes of the tensor */
	mat3x3 basis()
	{
		return _axisMatrix;
	}

	mat3x3 getTensor()
	{
		return _tensor;
	}

	double isotropicAverage()
	{
		return _isotropicAverage;
	}

	double anisotropyExtent();
	double smallness();
	vec3 longestAxis();
private:
	std::vector<vec3> _points;
	mat3x3 _tensor;
	void findPrincipleAxes();

	double _isotropicAverage;
	mat3x3 _axisMatrix;
	vec3 _axes[3];
};

#endif /* Anisotropicator_hpp */
