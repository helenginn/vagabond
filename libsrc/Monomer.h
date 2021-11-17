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
#include <hcsrc/FileReader.h>

/**
 * \class Monomer
 * \brief A class containing references to the atoms contained within a
 * residue of a Polymer.
 */

class Monomer : public AtomGroup
{
public:
	Monomer();
	virtual ~Monomer() {};
	void setup();
	void tieAtomsUp();

	bool isAfterAnchor();
	virtual void removeAtom(AtomPtr atom);

	MonomerPtr shared_from_this()
	{
		AtomGroupPtr group = AtomGroup::shared_from_this();
		return boost::static_pointer_cast<Monomer>(group);
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
	
	void empty();

	void setPolymer(PolymerPtr poly);

	PolymerPtr getPolymer()
	{
		return _myPolymer.lock();
	}

	virtual void addAtom(AtomPtr atom);

	virtual void setStream(std::ostream *str);

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
	
	AtomList getPhiAtoms();
	AtomList getPsiAtoms();

	std::string getResCode();

	virtual std::string getClassName()
	{
		return "Monomer";
	}

	virtual std::string getParserIdentifier()
	{
		return _identifier + i_to_str(_residueNum);
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void postParseTidy();
	virtual void linkReference(BaseParserPtr object, std::string category);
	
	static double vsRefine(void *object);
	void refine(CrystalPtr target,
	            RefinementType rType);

private:
	std::string _identifier; // e.g. three-letter code
	int _residueNum; // number in protein sequence including missing ones.

	PolymerWkr _myPolymer;
	std::vector<ModelPtr> _models;
	BackbonePtr _backbone;
	SidechainPtr _sidechain;

};

#endif /* defined(__vagabond__Monomer__) */
