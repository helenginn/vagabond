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

class Monomer : public std::enable_shared_from_this<Monomer>
{
public:
	Monomer();
	void setup();
	void tieAtomsUp();

	void setIdentifier(std::string idString)
	{
		_identifier = idString;
	}

	std::string getIdentifier()
	{
		return _identifier;
	}

	void setResidueNum(int n)
	{
		_residueNum = n;
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

	void addAtom(AtomPtr atom);

	void addModel(ModelPtr model)
	{
		_models.push_back(model);
	}

private:
	std::string _identifier; // e.g. three-letter code
	int _residueNum; // number in protein sequence including missing ones.

	PolymerWkr _myPolymer;
	std::vector<AtomPtr> _atoms;
	std::vector<ModelPtr> _models;
	BackbonePtr _backbone;
	SidechainPtr _sidechain;

};

#endif /* defined(__vagabond__Monomer__) */
