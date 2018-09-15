//
//  ExplicitModel.h
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__ExplicitModel__
#define __vagabond__ExplicitModel__

#include "Model.h"

/**
 * \class ExplicitModel
* \brief Abstract class providing a template for models which rely on
* explicit atom position calculations.
*
* These will include the Bond models and the anchor position for a polymer.
 */

class ExplicitModel : public Model
{
public:
	virtual ~ExplicitModel() {};

	void addRealSpacePositions(FFTPtr real, vec3 offset);

	virtual bool hasExplicitPositions()
	{
		return true;
	}
	FFTPtr makeDistribution();
protected:
	FFTPtr makeRealSpaceDistribution();

private:
};

#endif
