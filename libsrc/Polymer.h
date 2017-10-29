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

class Polymer :
public Molecule
{
public:
	Polymer()
	{
		_dampening = 0.05;
		_sideDampening = 0.05;
		_sideKick = 0;
		_anchorNum = 0;
		_totalMonomers = 0;
	}

	void closenessSummary();
	void addMonomer(MonomerPtr monomer);
	virtual void summary();
	virtual void tieAtomsUp();
	virtual void refine(CrystalPtr target, RefinementType rType);
	virtual void makePDB(std::string filename, PDBType pdbType, CrystalPtr crystal);
	virtual void graph(std::string graphName);
	virtual void differenceGraphs(std::string graphName, CrystalPtr diffCryst);

	static double getBackboneDampening(void *object);
	static void setBackboneDampening(void *object, double value);

	static double getSidechainDampening(void *object);
	static void setSidechainDampening(void *object, double value);

	static void setInitialKick(void *object, double value);
	static double getInitialKick(void *object);

	static double getSideKick(void *object);
	static void setSideKick(void *object, double value);

	void scaleFlexibilityToBFactor(CrystalPtr target);
	void scaleSidechainsToBFactor();
	void minimiseCentroids();
	void minimiseRotations();
	void downWeightResidues(int start, int end, double value);

	ModelPtr getAnchorModel();
	void changeAnchor(int num);
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

	long monomerCount()
	{
//		return _monomers.size();
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
	double _dampening;
	double _sideDampening;
	double _sideKick;
	double _totalMonomers;


};

#endif /* defined(__vagabond__Polymer__) */
