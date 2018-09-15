//
//  Anchor.h
//  vagabond
//
//  Created by Helen Ginn 2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __vagabond__Anchor__
#define __vagabond__Anchor__

#include "ExplicitModel.h"

/**
 * \class Anchor
 * \brief Model for an anchor point in a tied-up polymer. */

class Anchor : public ExplicitModel
{
public:
	Anchor(AbsolutePtr absolute);
	virtual ~Anchor() {};

	std::vector<BondSample> *getManyPositions();

	/** Returns the offsets for an anchor residue on which a molecule may
	* calculate translations and offsets.
	* \return x, y, z values centred around the origin.
	* */
	std::vector<vec3> getSphereAngles()
	{
		return _sphereAngles;
	}

	void setOccupancies(std::vector<double> occ)
	{
		_occupancies = occ;
	}
	
	double getBFactor()
	{
		return _bFactor;
	}
	
	void setBFactor(double b)
	{
		_bFactor = b;
	}
	
	AtomPtr getAtom()
	{
		return _atom.lock();
	}
protected:
	virtual std::string getClassName()
	{
		return "Anchor";
	}
	
private:
	vec3 _position;
	double _bFactor;
	AtomWkr _atom;

	std::vector<vec3> _sphereAngles;
	std::vector<double> _occupancies;
};

#endif
