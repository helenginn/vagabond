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
	Anchor();
	virtual ~Anchor() {};

	std::vector<BondSample> *getManyPositions(void *object = NULL);

	/** Returns the offsets for an anchor residue on which a molecule may
	* calculate translations and offsets.
	* \return x, y, z values centred around the origin.
	* */
	std::vector<vec3> getSphereAngles()
	{
		return _sphereAngles;
	}

	void setNeighbouringAtoms(AtomPtr nAtom, AtomPtr cAtom);

	AtomPtr getOtherAtom(AtomPtr calling);

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
	
	virtual std::string getParserIdentifier();

	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category);
private:
	void createStartPositions(Atom *callAtom);

	double _bFactor;
	AtomWkr _atom;
	AtomWkr _nAtom, _cAtom;
	
	vec3 _nDir;
	vec3 _cDir;

	std::vector<vec3> _sphereAngles;
	std::vector<double> _occupancies;
};

#endif