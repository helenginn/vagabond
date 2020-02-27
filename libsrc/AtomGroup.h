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

#define BUFFER_REGION (2.0)

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
	AtomList findAtomByNum(std::string atomType, int atomNum);

	double scoreWithMap(ScoreType scoreType, CrystalPtr crystal, 
	                    std::string plot = "", unsigned int flags = 0);

	
	static double scoreWithMapGeneral(MapScoreWorkspace *workspace,
	                                  bool plot = false);


	void makeBackboneTwists(ExplicitModelPtr applied);

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
	void addAtomsFrom(std::vector<AtomPtr> group);

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
	vec3 initialCentroid();

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

	size_t totalElements();
	
	ElementPtr element(size_t which)
	{
		return _elements[which];
	}
	
	void addToMap(VagFFTPtr fft, mat3x3 real2frac);

	/** Adds atoms to a map where the voxel morphology is cubic, with a
	  * given offset specified which is subtracted from each atom
	  * position. */
	void addToCubicMap(VagFFTPtr scratchFull);

	/** Prepares a cubic map to add the AtomGroup to, including adjustment
	 *  of the offset to place midpoint of the group of atoms at the midpoint 
	 * 	of the map. Returns the minimum Atom position as well.
	 *  \return offset which should be applied to each atom when calling
	 *  AtomGroup::addToCubicMap */
	void prepareCubicMap(VagFFTPtr *scratchFull, vec3 min, vec3 max, 
	                     bool cc = false);

	void setTargetRefinement(CrystalPtr target, RefinementType rType);
	virtual void refine(CrystalPtr target, RefinementType rType);
	void setWeighting(double value);
	int conformerCount();
	std::string conformer(size_t i);
	int conformer(std::string conf);
	
	/** Instructs the models of all the atoms inside to propagate a change of
	* parameters. See also: Model::propagateChange(). */
	void propagateChange();
	void propagateChangeExceptAnchor();
	void refreshPositions(bool quick = true);

	AtomPtr getClosestAtom(CrystalPtr crystal, vec3 pos);
	std::vector<AtomPtr> getHydrogenBonders();

	std::vector<AtomPtr> getAtoms()
	{
		return _atoms;
	}
	
	static double recalculatePositions(void *obj);
	
	void refreshBondAngles();
	virtual AtomList topLevelAtoms();
	
	void sort();
	void saveAtomPositions();

	void boundingMonomers(int *begin, int *end);
	
	void setName(std::string name)
	{
		_name = name;
	}
	
	void addToName(std::string name)
	{
		_name += name;
	}
	
	std::string getName()
	{
		return _name;
	}

	bool isTied()
	{
		return _beenTied;
	}

protected:
	virtual bool shouldRefineAtom(AtomPtr atom) { return true; };
	int _timesRefined;

	virtual std::string getClassName()
	{
		return "AtomGroup";
	}

	virtual std::string getParserIdentifier()
	{
		return "AtomGroup_" + _name;
	}
	
	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void linkReference(BaseParserPtr object, std::string category);
	std::vector<AtomPtr> _atoms;
	std::string _name;
private:
	static void plotCoordVals(std::vector<CoordVal> &vals, bool difference,
	                          double cutoff, std::string filename);
	void xyzLimits(vec3 *min, vec3 *max);

	static double scoreFinalMap(MapScoreWorkspace *workspace, bool plot,
	                            bool first);

	static double scoreFinalValues(std::vector<double> &xs,
	                               std::vector<double> &ys,
                                   std::vector<double> &weights,
	                               ScoreType scoreType,
                                   unsigned int flags);

	MonomerWkr _monomer;

	bool _beenTied;
	CrystalPtr _target;
	RefinementType _rType;

	void privateRefine(); 
	std::map<std::string, size_t> conformerMap();

	std::vector<ElementPtr> _elements;
	std::map<ElementPtr, FFTPtr> _eleScratch;
	size_t _scratchDims[3];
};

#endif /* defined(__vagabond__AtomGroup__) */
