//
//  Molecule.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Molecule__
#define __vagabond__Molecule__

#include <stdio.h>
#include <vector>
#include "mat3x3.h"
#include "shared_ptrs.h"
#include "Sampler.h"
#include "AtomGroup.h"
#include <string>
#include <iostream>

struct vec3;

/**
 * \class Molecule
 * \brief The highest level group of atoms directly referenced by a Crystal object.
 */

class Molecule : public AtomGroup
{
public:
	Molecule();

	void addToMap(FFTPtr fft, mat3x3 _real2frac, bool mask = false);

	virtual void summary();
	virtual void tieAtomsUp();
	virtual void refine(CrystalPtr target, RefinementType rType);
	void tiedUpScattering(double *tied, double *all);
	virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal);
	virtual void graph(std::string) {};
	virtual void differenceGraphs(std::string, CrystalPtr) {};
	void resetInitialPositions();
	void setAnchors();
	void makePowderList();

	virtual void reportParameters();

	double getAbsoluteBFacSubtract()
	{
		return _absoluteBFacSubtract;
	}

	void setAbsoluteBFacSubtract(double subtract)
	{
		_absoluteBFacSubtract = subtract;
	}

	double getAbsoluteBFacMult()
	{
		return _absoluteBFacMult;
	}

	void setAbsoluteBFacMult(double mult);

	void setChainID(std::string chain)
	{
		_chainID = chain;
	}

	std::string getChainID()
	{
		return _chainID;
	}

	std::vector<vec3> getCentroidOffsets()
	{
		return _centroidOffsets;
	}

	std::vector<mat3x3> getRotationCorrections()
	{
		return _rotations;
	}

	std::vector<vec3> getRotationCentres()
	{
		return _centroids;
	}

	std::vector<vec3> getTransTensorOffsets()
	{
		return _transTensorOffsets;
	}

	std::vector<mat3x3> getExtraRotations();

	vec3 getExtraRotationCentre()
	{
		return _rotationCentre;
	}

	virtual std::string getClassName()
	{
		return "Molecule";
	}

	bool isPolymer()
	{
		return (getClassName() == "Polymer");
	}

	MoleculePtr shared_from_this()
	{
		AtomGroupPtr groupPtr = AtomGroup::shared_from_this();
		return ToMoleculePtr(groupPtr);
	}

	static void setRotCentreZ(void *object, double value)
	{
		static_cast<Molecule *>(object)->_rotationCentre.z = value;
	}

	static double getRotCentreZ(void *object)
	{
		return static_cast<Molecule *>(object)->_rotationCentre.z;
	}

	static void setRotCentreY(void *object, double value)
	{
		static_cast<Molecule *>(object)->_rotationCentre.y = value;
	}

	static double getRotCentreY(void *object)
	{
		return static_cast<Molecule *>(object)->_rotationCentre.y;
	}

	static void setRotCentreX(void *object, double value)
	{
		static_cast<Molecule *>(object)->_rotationCentre.x = value;
	}

	static double getRotCentreX(void *object)
	{
		return static_cast<Molecule *>(object)->_rotationCentre.x;
	}

	std::vector<AtomPtr> getCloseAtoms(AtomPtr one, double tol, bool cache = false);
	
	void clearCloseCache()
	{
		_closeishAtoms.clear();
	}

protected:
	virtual std::string getParserIdentifier()
	{
		return "chain_" + _chainID; 
	}

	virtual void addProperties();
	virtual void postParseTidy();    

	std::vector<vec3> _centroidOffsets;
	std::vector<vec3> _centroids; // after offset correction
	std::vector<mat3x3> _rotations;

	std::vector<vec3> _transTensorOffsets;
	std::vector<mat3x3> _extraRotationMats; // currently unused

	virtual void calculateExtraRotations() {};

	// this axis calculates the angular response to the reaction sphere
	vec3 _magicRotAxis;

	// this axis is that of the rotation matrices applied to the structure
	vec3 _rotationAxis;
	vec3 _rotationCentre;

	double _rotationAngle;

	void setChangedRotation()
	{
		_changedRotations = true;
	}

private:
	double _absoluteBFacSubtract;
	double _absoluteBFacMult;

	bool _changedRotations;
	std::string _chainID;
	std::vector<AtomPtr> _closeishAtoms;
};

#endif /* defined(__vagabond__Molecule__) */
