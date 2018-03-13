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
#include "Parser.h"

class Anisotropicator;

class Absolute : public Model
{
public:
	Absolute(vec3 pos, double bFac, std::string element, double occValue);
	Absolute();

	// Model virtual functions:
	virtual std::vector<BondSample> *getManyPositions();
	virtual FFTPtr getDistribution(bool = false, int new_n = -1);
	virtual vec3 getAbsolutePosition()
	{
		return _position;
	}

	virtual void addToMolecule(MoleculePtr molecule);
	virtual void addToMonomer(MonomerPtr monomer);
	virtual mat3x3 getRealSpaceTensor();
	virtual void getAnisotropy(bool withKabsch);

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

	virtual void propagateChange(int depth = -1)
	{
		_calculated = false;
	}

	void setAlternativeConformerName(std::string conformer)
	{
		_conformer = conformer;
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
		_tensor = make_mat3x3();
		_tensor.vals[0] = bfac;
		_tensor.vals[4] = bfac;
		_tensor.vals[8] = bfac;
	}

	void setTensor(mat3x3 tensor, CrystalPtr crystal);

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

	virtual double getMeanSquareDeviation();

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

	std::vector<vec3> getSphereAngles()
	{
		return _sphereAngles;
	}

	void setAnchorPoint()
	{
		_isOfManyPositions = true;
	}
protected:
	static double getExpValue(void *object, double x, double y, double z);

	virtual std::string getParserIdentifier()
	{
		return "absolute_" + i_to_str(getAtom()->getAtomNum());
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category) {};
	virtual void linkReference(ParserPtr object, std::string category);
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
	std::string _conformer;
	std::vector<BondSample> _bondSamples;
	std::vector<vec3> _sphereAngles;

	vec3 _position;
	double _bFactor;
	bool _isOfManyPositions;

	void makeAtom();
};

#endif /* defined(__vagabond__Absolute__) */
