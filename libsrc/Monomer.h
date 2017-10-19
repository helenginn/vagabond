//
//  Monomer.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Monomer__
#define __vagabond__Monomer__

#include <stdio.h>
#include <string>
#include "shared_ptrs.h"
#include <vector>
#include "Sidechain.h"
#include "AtomGroup.h"

class Monomer : public AtomGroup
{
public:
	Monomer();
	void setup();
	void tieAtomsUp();

	/* For global modification of entire structure */
	void setBackboneDampening(double value);
	void setSidechainDampening(double value);

	bool isAfterAnchor();

	MonomerPtr shared_from_this()
	{
		AtomGroupPtr group = AtomGroup::shared_from_this();
		return std::static_pointer_cast<Monomer>(group);
	}

	void setIdentifier(std::string idString)
	{
		_identifier = idString;
	}

	std::string getIdentifier()
	{
		return _identifier;
	}

	BackbonePtr getBackbone()
	{
		return _backbone;
	}

	SidechainPtr getSidechain()
	{
		return _sidechain;
	}

	void setResidueNum(int n)
	{
		_residueNum = n;
		_sidechain->setResNum(_residueNum);
	}

	int getResidueNum()
	{
		return _residueNum;
	}

	void setPolymer(PolymerPtr poly)
	{
		_myPolymer = poly;
	}

	PolymerPtr getPolymer()
	{
		return _myPolymer.lock();
	}

	virtual void addAtom(AtomPtr atom);

	void addModel(ModelPtr model)
	{
		_models.push_back(model);
	}

	int modelCount()
	{
		return _models.size();
	}

	ModelPtr model(int i)
	{
		return _models[i];
	}

	void setKick(double value, bool beforeAnchor);
	double getKick();

	void setSideKick(double value);
private:
	std::string _identifier; // e.g. three-letter code
	int _residueNum; // number in protein sequence including missing ones.

	PolymerWkr _myPolymer;
	std::vector<ModelPtr> _models;
	BackbonePtr _backbone;
	SidechainPtr _sidechain;

};

#endif /* defined(__vagabond__Monomer__) */
