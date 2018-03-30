//
//  Model.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Model__
#define __vagabond__Model__

#include "shared_ptrs.h"
#include <stdio.h>
#include <vector>
#include "Distributor.h"
#include <mutex>
#include "Parser.h"

typedef struct
{
	mat3x3 basis;   /* Defines bond axis of previous bond */
	vec3 start;     /* position of last minor */
	vec3 old_start; /* position of torsion-defining atom */
	double torsion; /* Defines torsion of next atom */
	double occupancy;
} BondSample;

// Anything which is capable of predicting electron positions.
//

class Model : public boost::enable_shared_from_this<Model>, public Distributor, public Parser
{
public:
	Model();

	virtual FFTPtr getDistribution(bool quick = false, int n = -1) = 0;

	virtual void addToMonomer(MonomerPtr monomer);
	virtual void addToMolecule(MoleculePtr) {};

	virtual std::string getClassName() = 0;

	/** Actual mean position of blurred positions including whole-molecule
	* 	deviations. */
	virtual vec3 getAbsolutePosition()
	{
		return _absolute;
	}

	/** Mean position of blurred positions without including whole-molecule
	* 	deviations. Should not be used for final atom position calculations. */
	virtual std::vector<BondSample> *getManyPositions() = 0;

	/** Individual positions including whole-molecule deviations. */
	std::vector<vec3> polymerCorrectedPositions();

	/** Positions and associated data including whole-molecule deviations.
	* 	Will return from cache if not flagged to recalculate. */
	virtual std::vector<BondSample> getFinalPositions();

	virtual double getEffectiveOccupancy() { return 1; }

	virtual double getMeanSquareDeviation() = 0;
	
	/**
	* Get the tensor as described by Vagabond model. May be different to
	* that rederived from the Vagabond models. 
	* \return symmetric matrix with 6 unique values
	*/
	virtual mat3x3 getRealSpaceTensor();
	
	/**
	* Reports the standard deviation of the longest dimension of the 
	* anisotropic tensor corresponding to this distribution.
	* \return standard deviation (Angstroms)
	*/
	double biggestStdevDim();
	int fftGridLength();

	bool hasMolecule()
	{
		return !_molecule.expired();
	}

	MoleculePtr getMolecule()
	{
		return _molecule.lock();
	}

	void setMolecule(MoleculePtr mole)
	{
		_molecule = mole;
	}

	/**
	* For recursive models (e.g. Bonds) downstream models are flagged
	* to indicate that they need to be recalculated before using stored
	* position information.
	* \param maximum depth to flag propagation in number of models.
	* \param refresh immediately recalculate atom positions if true.
	*/
	virtual void propagateChange(int depth = -1, bool refresh = false);

	bool isBond()
	{
		return (getClassName() == "Bond");
	}

	bool isAbsolute()
	{
		return (getClassName() == "Absolute");
	}

	bool isAnchor()
	{
		return (getClassName() == "Anchor");
	}
	
	void setPolymerChanged()
	{
		_recalcFinal = true;	
	}

	vec3 longestAxis();
	std::vector<vec3> fishPositions();
protected:
	mat3x3 _realSpaceTensor;

	/** Molecule which provides whole-molecule offsets/rotations */
	MoleculeWkr _molecule;

	/** Caching of an atom's position. */
	vec3 _absolute;

	/** A record of the final positions just for GUI display */
	std::vector<vec3> _finalPositions;

	/** 
	* 	If using mutex, lock while modifying positions to prevent interference
	* from GUI */
	static bool _useMutex;

	vec3 _longest;
	double _anisotropyExtent;

	virtual void getAnisotropy(bool withKabsch);
	double anisotropyExtent(bool withKabsch = false);
	double _isotropicAverage;
	bool _recalcFinal;
	std::vector<BondSample> _finalSamples;

	virtual std::string getParserIdentifier()
	{
		return "model"; 
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr, std::string) {};
private:
	std::mutex guiLock;
};

#endif /* defined(__vagabond__Model__) */
