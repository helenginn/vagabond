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
#include "Bond.h"

class Absolute : public Model
{
public:
	Absolute(vec3 pos, double bFac, std::string element, double occValue);
	Absolute();

// Model virtual functions:
	virtual std::vector<BondSample> *getManyPositions(BondSampleStyle style);
	virtual FFTPtr getDistribution(bool quick = false);
	virtual vec3 getStaticPosition();
	virtual vec3 getAbsolutePosition()
	{
		return getStaticPosition();
	}
	
	virtual void addToMolecule(MoleculePtr molecule);
	virtual void addToMonomer(MonomerPtr monomer);

	void setIdentity(int resNumValue, std::string chainID,
					 std::string resName, std::string atomName, int atomNum)
	{
		_resNum = resNumValue;
		_chainID = chainID;
		_resName = resName;
		_atomName = atomName;
		_atomNum = atomNum;

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

	void setBFactor(double bfac)
	{
		_bFactor = bfac;
	}

	void setTensor(mat3x3 tensor)
	{
		_tensor = tensor;
		_usingTensor = true;
		double x = sqrt(tensor.vals[0]);
		double y = sqrt(tensor.vals[4]);
		double z = sqrt(tensor.vals[8]);
		getAtom()->setInitialAnisoBs(x, y, z);
	}

	double getBFactor()
	{
		return _bFactor;
	}

	void setHeteroAtom(bool hetatm)
	{
		_hetatm = hetatm;
	}

	bool isHeteroAtom()
	{
		return _hetatm;
	}

	static void setPosX(void *object, double x)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_position.x = x;
	}

	static void setPosY(void *object, double y)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_position.y = y;
	}

	static double getPosZ(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_position.z;
	}

	static double getPosX(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_position.x;
	}

	static double getPosY(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_position.y;
	}

	static void setPosZ(void *object, double z)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_position.z = z;
	}

	virtual double getMeanSquareDeviation(double target = -1, int index = -1);

	static double getB(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_bFactor;
	}

	static void setB(void *object, double b)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_bFactor = b;
		abs->_calculated = false;
	}

	virtual std::string getClassName()
	{
		return "Absolute";
	}

	AtomPtr getAtom()
	{
		return _atom;
	}

	void addNextAtom(AtomPtr atom)
	{
		_nextAtoms.push_back(atom);
	}

	AtomPtr getNextAtom(int i)
	{
		return _nextAtoms[i].lock();
	}

	long nextAtomCount()
	{
		return _nextAtoms.size();
	}
protected:
	static double getExpValue(void *object, double x, double y, double z);

private:
	AtomPtr _atom;
	std::vector<AtomWkr> _nextAtoms;
	std::string _element;
	double _occupancy;
	std::string _chainID, _resName, _atomName;
	int _resNum, _atomNum;
	bool _hetatm;
	bool _usingTensor;
	mat3x3 _tensor;
	std::vector<BondSample> _bondSamples;

	vec3 _position;
	double _bFactor;

	void makeAtom();
};

#endif /* defined(__vagabond__Absolute__) */
