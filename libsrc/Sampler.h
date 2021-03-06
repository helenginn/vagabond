//
//  Sampler.h
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Sampler__
#define __vagabond__Sampler__

#include <stdio.h>
#include <vector>
#include "MapScoreWorkspace.h"
#include "mat3x3.h"
#include "shared_ptrs.h"
#include "RefinementStrategy.h"
#include <string>
#include <map>

/**
 * \class Sampler
 * \brief Organises a target function and samples parameters to optimise
 * that target function.
 */

/** Flags to set refinement strategies for a protein chain. */
typedef enum
{
	RefinementCrude = 0, /** Refinement from far soln against ED */
	RefinementFine = 1, /** Refinement against electron density */
	RefinementSidechain = 2, /** Sidechains against electron density */
	RefinementModelPos = 4, /** Positions to PDB positions */
	RefinementCentroid = 6, /** Refine centroid to PDB centroid */
	RefinementWaterNetwork = 7,
	RefinementMouse = 8, /** Refinement to mouse coordinates */
	RefinementSavedPos = 9, /** Refinement to previously saved positions */
	RefinementSidePos = 10, /** Only position params to sidechain density */
} RefinementType; 


typedef std::map<ParamOptionType, double> ParamMap;

class SVDBond;

class Sampler
{
public:
	Sampler();
	
	virtual ~Sampler()
	{

	}

	/**
	* 	 Add an atom from the sensitive area in real space
	*/
	void addSampled(AtomPtr atom);
	
	/**
	*  Add a bunch of atoms from the sensitive area in real space	
	*/
	void addSampled(std::vector<AtomPtr> atoms);

	/**
	* Add torsion angle to sampled parameters
	* \param bond bond to vary torsion around
	* \param range step size for angle in degrees
	* \param interval unused usually, but could be stopping point. 	 
	*/
	void addTorsion(BondPtr bond, double range, double interval);

	void addTwist(BondPtr bond, double range, double interval);

	/**
	* Add torsion blur (kick function) to sampled parameters
	* \param bond bond to vary kick for
	* \param range step size in degree offset per degree
	* \param interval unused usually, but could be stopping point. 	 
	*/
	void addKick(BondPtr bond, double range, double interval);
	void addBondLength(BondPtr bond, double range, double interval);
	void addBendAngle(BondPtr bond, double range, double interval,
	                  bool circlePortionOnly = false);
	void addOccupancy(BondPtr bond, double range, double interval);

	/** Add sampled atoms from a given atom group and conformer name */
	void addSampledAtoms(AtomGroupPtr group, std::string conformer = "");
	
	/** Add anchor position for an atom controlled by Anchor class */
	void addAnchorParams(AnchorPtr anch);

	/** Add absolute B factor for an atom controlled by Absolute class */
	void addAbsoluteBFactor(AbsolutePtr abs, double range, double interval);
	void addMagicAngle(BondPtr bond, double range, double interval);
	void setCrystal(CrystalPtr crystal);
	bool sample(bool clear = true);

	void setupGrid();
	void setupNelderMead();
	void setupStepSearch();

	void reportInDegrees()
	{
		_strategy->reportInDegrees();
	}

	void setCycles(int cycles)
	{
		_cycles = cycles;
	}

	static double score(void *object)
	{
		return static_cast<Sampler *>(object)->getScore();
	}

	void setJobName(std::string newJob)
	{
		_jobName = newJob;
	}

	std::string getJobName()
	{
		return _jobName;
	}
	
	void convertToSVD(bool conv)
	{
		_convert = conv;
	}

	void setSilent(bool silent = true)
	{
		_strategy->setSilent(silent);
		_silent = silent;
	}

	void setVerbose(bool verbose = true)
	{
		_verbose = true;
	}
	
	int sampleSize()
	{
		return _sampled.size();
	}

	/**
	* 	If mock is set, then the refinement is performed only to see the
	* 	change in target function and then the original values for the
	* 	parameters are restored.
	*/
	void setMock()
	{
		_mock = true;
	}
	
	bool didChange()
	{
		return _changed;
	}
	
	void saveScore();
	
	double getImprovement()
	{
		return _improv;
	}

	void setScoreType(ScoreType type)
	{
		_scoreType = type;
	}

	void clearParams()
	{
		_params.clear();
	}

	void addParamType(ParamOptionType type, double value)
	{
		_params[type] = value; 
	}

	size_t paramCount()
	{
		return _params.size();
	}

	void copyParams(SamplerPtr sampler)
	{
		sampler->_params = _params;
	}

	int hasParameter(ParamOptionType type);
	void checkOccupancyAndAdd(BondPtr bond);

	void addCustomParameter(void *object, Getter getter, Setter setter,
                                 double range, double interval,
                                 std::string name);


	std::vector<AtomGroupPtr> includingInRefinement()
	{
		return _includeForRefine;
	}


	void addIncludeForRefinement(AtomGroupPtr group)
	{
		_includeForRefine.push_back(group);
	}
protected:
	/** Create a torsion set which adds the primary chain (the backbone)

	 * 	appropriate bond parameters and the sensitive atoms. */	
	BondPtr setupThoroughSet(BondPtr bond, bool addBranches = true);
	mat3x3 _real2Frac;
	bool _excludeO;

	void addTwistShift(ExplicitModelPtr eModel, AtomGroupPtr clearGroup);
	int _refinedMagicAxisCount;
	virtual double getScore();
	void setupCloseAtoms();
	void setupScoreWithMap();
	void addAtomsForBond(BondPtr bond);

	double getParameter(ParamOptionType type)
	{
		return _params[type];
	}
	
	ParamMap getParams()
	{
		return _params;
	}
	
	RefinementStrategyPtr getStrategy()
	{
		return _strategy;
	}

	void clearIncludeForRefinements()
	{
		_includeForRefine.clear();
	}

	std::vector<AtomGroupPtr> _includeForRefine;
private:
	double constraint();
	void addParamsForBond(BondPtr bond, bool even = true);
	CrystalPtr _crystal;
	std::vector<BalancePtr> _balances;

	std::vector<AtomPtr> _sampled;
	std::vector<BondPtr> _bonds;
	ParamMap _params;

	bool _mock;
	bool _silent;
	bool _changed;
	bool _convert;
	double _improv;
	double _begin;
	bool _shouldSave;
	bool _verbose;
	int _cycles;

	SVDBond *_svd;

	std::string _jobName;
	ScoreType _scoreType;

	RefinementStrategyPtr _strategy;
	MapScoreWorkspace _workspace;
};


#endif /* defined(__vagabond__Sampler__) */
