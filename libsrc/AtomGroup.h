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
#include "MapScoreWorkspace.h"
#include "Sampler.h"
#include "Parser.h"
#include <map>
#include <climits>

/**
 * \class AtomGroup
 * \brief AtomGroup looks after the concept of any sensible group of Atoms,
 * and is often (but does not need to be) subclassed.
 */

class AtomGroup : public Sampler, public Parser
{
public:
	AtomGroupPtr shared_from_this()
	{
		return ToAtomGroupPtr(Parser::shared_from_this());
	}

	AtomGroup();
	virtual ~AtomGroup() {};
	AtomPtr findAtom(std::string atomType);
	AtomPtr findAtom(std::string atomType, std::string confID);
	AtomList findAtoms(std::string atomType);
	AtomGroupPtr subGroupForConf(int conf);
	AtomList findAtoms(std::string atomType, int resNum);

	double scoreWithMap(ScoreType scoreType, CrystalPtr crystal, 
	                    bool plot = false, unsigned int flags = 0);

	
	static double scoreWithMapGeneral(MapScoreWorkspace *workspace,
	                                  bool plot = false);


	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}

	MonomerPtr getMonomer()
	{
		return _monomer.lock();
	}

	virtual void addAtom(AtomPtr atom);
	
	virtual void removeAtom(AtomPtr atom);

	void addAtomsFrom(AtomGroupPtr group);

	AtomList beyondGroupAtoms(bool just_bottom = false);
	
	size_t atomCount()
	{
		return _atoms.size();
	}

	bool hasAtom(AtomPtr anAtom);

	AtomPtr atom(int i)
	{
		return _atoms[i];
	}
	
	vec3 centroid();

	double totalElectrons();
	double getAverageBFactor(bool initial = false);
	
	/** Average offset from initial PDB position */
	double getAverageDisplacement();

	std::string getPDBContribution(PDBType pdbType,
	                               CrystalPtr crystal = CrystalPtr(),
	                               int conformer = -1);

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
	std::string conformer(size_t i);
	int conformer(std::string conf);
	
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

	AtomPtr getClosestAtom(CrystalPtr crystal, vec3 pos);
	std::vector<AtomPtr> getHydrogenBonders();

	std::vector<AtomPtr> getAtoms()
	{
		return _atoms;
	}
	
	void refreshBondAngles();
	virtual AtomList topLevelAtoms();
	
	void saveAtomPositions();
	
	std::vector<AtomGroupPtr> includingInRefinement()
	{
		return _includeForRefine;
	}
protected:
	virtual bool shouldRefineAtom(AtomPtr atom) { return true; };
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
	std::vector<AtomPtr> _atoms;
	std::vector<AtomGroupPtr> _includeForRefine;
private:
	static void plotCoordVals(std::vector<CoordVal> &vals, bool difference,
	                          double cutoff, std::string filename);
	static FFTPtr prepareMapSegment(CrystalPtr crystal,
	                                std::vector<AtomPtr> selected,
	                                mat3x3 *basis, vec3 *ave);

	static double scoreFinalMap(CrystalPtr crystal, FFTPtr segment,
	                            bool plot, ScoreType scoreType,
	                            vec3 ave, unsigned int flags = 0);

	static double scoreFinalValues(std::vector<double> xs,
	                               std::vector<double> ys,
	                               ScoreType scoreType,
                                   unsigned int flags);

	MonomerWkr _monomer;

	bool _beenTied;
	CrystalPtr _target;
	RefinementType _rType;

	void privateRefine(); 
	std::map<std::string, size_t> conformerMap();

};

#endif /* defined(__vagabond__AtomGroup__) */
