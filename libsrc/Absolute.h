//
//  Absolute.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Absolute__
#define __vagabond__Absolute__

#include <stdio.h>
#include "Model.h"
#include "vec3.h"
#include <string>
#include "FileReader.h"
#include "Distributor.h"

class Absolute : public Model, public Distributor
{
public:
	Absolute(vec3 pos, double bFac, std::string element, double occValue);


// Model virtual functions:
	virtual FFTPtr getDistribution();
	virtual vec3 getPosition();
	virtual void addToMolecule(MoleculePtr molecule);
	virtual void addToMonomer(MonomerPtr monomer);

	void setIdentity(int resNumValue, std::string chainID,
					 std::string resName, std::string atomName)
	{
		_resNum = resNumValue;
		_chainID = chainID;
		_resName = resName;
		_atomName = atomName;

		trim(_chainID); trim(_atomName);

		trim(_resName); to_lower(_resName);
	}

	int getResNum()
	{
		return _resNum;
	}

	std::string getResName()
	{
		return _resName;
	}

	std::string getChainID()
	{
		return _chainID + (_hetatm ? "_hetatm" : "");
	}

	void setHeteroAtom(bool hetatm)
	{
		_hetatm = hetatm;
	}

	bool isHeteroAtom()
	{
		return _hetatm;
	}

	virtual std::string getClassName()
	{
		return "Absolute";
	}
protected:
	static double getExpValue(void *object, double x, double y, double z);

private:
	AtomPtr _atom;
	std::string _element;
	double _occupancy;
	std::string _chainID, _resName, _atomName;
	int _resNum;
	bool _hetatm;

	vec3 _position;
	double bFactor;

	void makeAtom();
};

#endif /* defined(__vagabond__Absolute__) */
