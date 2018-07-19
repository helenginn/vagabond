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

/**
 * \class Model
* \brief Abstract class providing a template for anything which is capable of
* determining an atom's position and distribution.
*
* The two major classes of Model (as of March 2018) are Absolute and Bond
* models. The Absolute class is akin to an entry in a PDB file, whereas the
* Bond class is for Vagabond-specific recursive refinement.
*
* The Model abstract class does look after a few things, like whole-molecule
* rotations and translations, caching of the absolute position(s) and allowing
* the GUI to fish out the positions for display.
 */

/** \struct BondSample
 *  \brief Transfers information between Bonds for calculation of positional
 *  information.
 */
typedef struct
{
	mat3x3 basis;     /**< Defines bond axis of previous bond */
	vec3 start;       /**< position of last minor */
	vec3 old_start;   /**< position of torsion-defining atom */
	double torsion;   /**< Defines torsion of next atom */
	double occupancy; /**< Relative occupancy (usually 1) */
	double kickMult;  /**< Prelim work - some sub-groups can kick more */
} BondSample;

class Model : public Distributor, public Parser
{
public:
	ModelPtr shared_from_this()
	{
		return ToModelPtr(Parser::shared_from_this());
	}

	Model();
	virtual ~Model() {};

	virtual FFTPtr makeDistribution() = 0;
	FFTPtr getDistribution();

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

	/** Occupancy for a given atom after all modifiers applied.
	* 	\return value between 0 and 1. */
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
	
	/**
	* Suggests a grid length to use for an FFT containing this atom.
	* \return integer grid length for x=y=z.
	*/
	int fftGridLength();
	
	virtual AtomPtr getAtom() = 0;

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
		_recalcDist = true;
	}
	
	void setRecalculateDist()
	{
		_recalcDist = true;
	}
	
	vec3 getSpecificPosition(int i)
	{
		if (i > _finalPositions.size()) return empty_vec3();
		return getFinalPositions()[i].start;
	}
	
	void addRealSpacePositions(FFTPtr real, vec3 offset);

	vec3 longestAxis();
	double smallness();
	std::vector<vec3> fishPositions();
	
	virtual bool hasExplicitPositions() = 0;
	void writePositionsToFile(std::string filename);
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
	double _smallness;
	FFTPtr _lastDistribution;

	virtual void getAnisotropy(bool withKabsch);
	double anisotropyExtent(bool withKabsch = false);
	double _isotropicAverage;
	bool _recalcFinal;
	bool _recalcDist;
	std::vector<BondSample> _finalSamples;
	
	FFTPtr makeRealSpaceDistribution();

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
