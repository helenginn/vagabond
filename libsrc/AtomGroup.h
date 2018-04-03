//
//  AtomGroup.h
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__AtomGroup__
#define __vagabond__AtomGroup__

#include <stdio.h>
#include <string>
#include <vector>
#include "shared_ptrs.h"
#include "Sampler.h"
#include "Parser.h"
#include <map>

/**
 * \class AtomGroup
 * \brief AtomGroup looks after the concept of any sensible group of Atoms,
 * and is often (but does not need to be) subclassed.
 */

class AtomGroup : public boost::enable_shared_from_this<AtomGroup>, public Sampler, public Parser
{
public:
	AtomGroup();
	AtomPtr findAtom(std::string atomType);
	AtomPtr findAtom(std::string atomType, std::string confID);
	AtomList findAtoms(std::string atomType);

	double scoreWithMap(ScoreType scoreType, CrystalPtr crystal, bool plot = false);
	static double scoreWithMapGeneral(ScoreType scoreType, CrystalPtr crystal,
	                                  bool plot = false,
	std::vector<AtomPtr> selected = std::vector<AtomPtr>());

	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}

	MonomerPtr getMonomer()
	{
		return _monomer.lock();
	}

	virtual void addAtom(AtomPtr atom)
	{
		std::vector<AtomPtr>::iterator it;
		it = std::find(_atoms.begin(), _atoms.end(), atom);

		if (it == _atoms.end())
		{
			_atoms.push_back(atom);
		}
	}
	
	void addAtomsFrom(AtomGroupPtr group);

	long atomCount()
	{
		return _atoms.size();
	}

	bool hasAtom(AtomPtr anAtom);

	AtomPtr atom(int i)
	{
		return _atoms[i];
	}

	double totalElectrons();
	double getAverageBFactor(bool initial = false);
	double getAverageDisplacement();

	std::string getPDBContribution(PDBType pdbType,
	                               CrystalPtr crystal = CrystalPtr());

	void setTied()
	{
		_beenTied = true;
	}

	virtual bool shouldRefineAngles()
	{
		return false;
	}

	int totalElectrons(int *fcWeighted);

	static double refine(void *object)
	{
		static_cast<AtomGroup *>(object)->privateRefine();
		return 0;
	}

	void setTargetRefinement(CrystalPtr target, RefinementType rType);
	virtual void refine(CrystalPtr target, RefinementType rType);
	void setWeighting(double value);
	int conformerCount();
	std::string conformer(int i);
	
	/** Instructs the models of all the atoms inside to propagate a change of
	* parameters. See also: Model::propagateChange(). */
	void propagateChange();
	void refreshPositions(bool quick = true);

	void clearIncludeForRefinements()
	{
		_includeForRefine.clear();
	}

	void addIncludeForRefinement(AtomGroupPtr group)
	{
		_includeForRefine.push_back(group);
	}
protected:
	virtual AtomList topLevelAtoms();
	int _timesRefined;

	bool isTied()
	{
		return _beenTied;
	}

	virtual std::string getClassName()
	{
		return "AtomGroup";
	}

	virtual std::string getParserIdentifier()
	{
		return "AtomGroupSomething";
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void linkReference(ParserPtr object, std::string category);
private:
	MonomerWkr _monomer;

	std::vector<AtomPtr> _atoms;

	bool _beenTied;
	CrystalPtr _target;
	RefinementType _rType;
	std::vector<AtomGroupPtr> _includeForRefine;

	void privateRefine(); 
	std::map<std::string, int> conformerMap();

};

#endif /* defined(__vagabond__AtomGroup__) */
