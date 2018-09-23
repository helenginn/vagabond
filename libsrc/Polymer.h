//
//  Polymer.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Polymer__
#define __vagabond__Polymer__

#include <stdio.h>
#include "shared_ptrs.h"
#include "Molecule.h"
#include <vector>
#include <map>
#include "Options.h"

/**
 * \class Polymer
 * \brief A subclass of Molecule which contains a series of Monomer objects,
 * forming a polymer chain.
 */

class FlexGlobal;

class Polymer : public Molecule
{
public:
	Polymer();
	virtual ~Polymer() {}

	void closenessSummary();
	void addMonomer(MonomerPtr monomer);
	virtual void summary();
	virtual void tieAtomsUp();
	void splitConformers();
	virtual void refine(CrystalPtr target, RefinementType rType);
	
	static void refineVScript(void *object, RefinementType rType);
	static double vsRefineSidechainsToDensity(void *object);
	static double vsRefinePositionsToPDB(void *object);
	static double vsRefineLocalFlexibility(void *object);
	static double vsRefineGlobalFlexibility(void *object);
	
	virtual void graph(std::string graphName);
	virtual void differenceGraphs(std::string graphName, CrystalPtr diffCryst);

	static void setInitialKick(void *object, double value);
	static double getInitialKick(void *object);

	static double vsSandbox(void *object);
	
	void scaleSidechainsToBFactor();
	void refineBackbone();
	void refineBackboneFrom(int position);
	static double vsRefineBackbone(void *object);
	static void vsRefineBackboneFrom(void *object, double position);
	static void vsMultiplyBackboneKick(void *object, double value);

	static void vsOmitResidues(void *object, double start, double end);
	static void vsUnomitResidues(void *object, double start, double end);
	
	void attachTargetToRefinement(RefinementStrategyPtr strategy,
	                              FlexGlobal &target, bool isotropy = false);

	virtual void reportParameters();
	void downWeightResidues(int start, int end, double value);

	void applyPolymerChanges();
	void refineToEnd(int monNum, CrystalPtr target, RefinementType rType);
	double refineRange(int start, int end, 
	                 CrystalPtr target, RefinementType rType);
	bool test();
	ExplicitModelPtr getAnchorModel();
	void findAnchorNearestCentroid();
	void hydrogenateContents();
	void checkChainContinuity();
	void setAnchor(int num)
	{
		_anchorNum = num;
	}

	int getAnchor()
	{
		return _anchorNum;
	}

	MonomerPtr getMonomer(int i)
	{
		if (_monomers.count(i))
		{
			return _monomers[i];
		}

		return MonomerPtr();
	}

	int monomerCount()
	{
		//        return _monomers.size();
		return _totalMonomers;
	}


	virtual std::string getClassName()
	{
		return "Polymer";
	}

	PolymerPtr shared_from_this()
	{
		return ToPolymerPtr(Molecule::shared_from_this());
	}

	void refineAnchorMovements();
	void refineLocalFlexibility();
	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void postParseTidy();
	
	AtomGroupPtr getAllBackbone();
	
	std::map<long, MonomerPtr>::iterator beginMonomer()
	{
		return _monomers.begin();
	}
	
	double getKickShift()
	{
		return _kickShift;
	}
	
	void setKickShift(double shift)
	{
		_kickShift = shift;
	}
protected:
	virtual double getScore()
	{
		propagateChange();
		return Sampler::getScore();
	}

private:
	void refineMonomer(MonomerPtr monomer, CrystalPtr target,
	                   RefinementType rType);

	std::map<long, MonomerPtr> _monomers;

	int _anchorNum;
	double _startB;
	double _kick;
	double _kickShift;
	int _totalMonomers;

	AtomGroupPtr _allBackbones;

};

#endif /* defined(__vagabond__Polymer__) */
