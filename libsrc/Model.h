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


class Model : public Distributor, public Parser
{
public:
	ModelPtr shared_from_this()
	{
		return ToModelPtr(Parser::shared_from_this());
	}

	Model();
	virtual ~Model() {};

	/** Return an atom probability distribution from this model. Needs
	 * 	to be reimplemented by non-abstract models. */
	virtual FFTPtr makeDistribution() = 0;
	FFTPtr getDistribution();

	virtual std::string getClassName() = 0;

	/** To return the mean square deviation as a B factor (i.e.,
	 * variance * 8 * M_PI^2. Needs to be reimplemented by non-abstract 
	 * models. */
	virtual double getMeanSquareDeviation() = 0;

	/** Actual mean position of blurred positions including whole-molecule
	* 	deviations. */
	virtual vec3 getAbsolutePosition()
	{
		return _absolute;
	}

	/** Occupancy for a given atom after all modifiers applied.
	* 	\return value between 0 and 1. */
	virtual double getEffectiveOccupancy() { return 1; }
	
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
	
	/** Return the directly controlled atom of this model. Needs
	 * reimplementing by non-abstract model classes. */
	virtual AtomPtr getAtom() = 0;

	/** Returns true if this model is associated with a molecule. */
	bool hasMolecule()
	{
		return !_molecule.expired();
	}

	/** Returns the molecule associated with this model. Will fail if
	 * molecule is not assigned. Check with hasMolecule() first. */
	MoleculePtr getMolecule()
	{
		return _molecule.lock();
	}

	void setMolecule(MoleculePtr mole)
	{
		_molecule = mole;
	}

	virtual void addToMonomer(MonomerPtr monomer);
	virtual void addToMolecule(MoleculePtr) {};

	/**
	* For recursive models (e.g. Bonds) downstream models are flagged
	* to indicate that they need to be recalculated before using stored
	* position information.
	* \param maximum depth to flag propagation in number of models.
	* \param refresh immediately recalculate atom positions if true.
	*/
	virtual void propagateChange(int depth = -1, bool refresh = false);
	
	virtual void refreshPositions() = 0;

	bool isBond()
	{
		return (getClassName() == "Bond");
	}

	bool isSponge()
	{
		return (getClassName() == "Sponge");
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
	
	virtual vec3 longestAxis();
	double smallness();
	
	virtual bool hasExplicitPositions() = 0;
	
	double getBFactor()
	{
		return getMeanSquareDeviation();
	}
	
	virtual std::string shortDesc() = 0;
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

	virtual void getAnisotropy(bool) {};
	virtual double anisotropyExtent(bool withKabsch = false);
	double _isotropicAverage;
	bool _recalcFinal;
	bool _recalcDist;
	
	virtual std::string getParserIdentifier()
	{
		return "model"; 
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr, std::string) {};
private:
};

#endif /* defined(__vagabond__Model__) */
