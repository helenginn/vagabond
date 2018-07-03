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

class Plucker;

/**
 * \class Molecule
 * \brief The highest level group of atoms directly referenced by a Crystal object.
 */

class Molecule : public AtomGroup
{
public:
	Molecule();
	virtual ~Molecule() {};

	void addToMap(FFTPtr fft, mat3x3 _real2frac, bool mask = false);

	virtual void addAtom(AtomPtr atom);
	
	virtual void summary();
	virtual void tieAtomsUp();
	virtual void refine(CrystalPtr target, RefinementType rType);
	void tiedUpScattering(double *tied, double *all);
	virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal,
	                            int conformer = -1);
	virtual void graph(std::string) {};
	virtual void differenceGraphs(std::string, CrystalPtr) {};
	void resetInitialPositions();
	void setAnchors();
	void makePowderList();
	void expandWaters();

	virtual void reportParameters();

	double getAbsoluteBFacSubtract()
	{
		return _absoluteBFacSubtract;
	}
	
	static double getAbsoluteBFacSubtract(void *object)
	{
		return static_cast<Molecule *>(object)->_absoluteBFacSubtract;
	}

	void setAbsoluteBFacSubtract(double subtract)
	{
		_absoluteBFacSubtract = subtract;
		std::cout << "Setting absolute B factor subtractor to " << subtract << std::endl;
		refreshBModels();
	}
	
	static void setAbsoluteBFacSubtract(void *object, double subtract);

	static void vsSetAbsoluteBFacSubtract(void *object, double value)
	{
		Parser *parser = static_cast<Parser *>(object);
		Molecule *molecule = dynamic_cast<Molecule *>(parser);
		molecule->setAbsoluteBFacSubtract(value);	
	}

	double getAbsoluteBFacMult()
	{
		return _absoluteBFacMult;
	}

	void refreshBModels();
	void setAbsoluteBFacMult(double mult);
	static void vsSetAbsoluteBFacMult(void *object, double mult)
	{
		Parser *parser = static_cast<Parser *>(object);
		Molecule *molecule = dynamic_cast<Molecule *>(parser);
		molecule->setAbsoluteBFacMult(mult);	
	}

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

	bool isWaterNetwork()
	{
		return (getClassName() == "WaterNetwork");
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
	std::vector<mat3x3> _extraRotationMats;

	virtual void calculateExtraRotations() {};

	// this axis calculates the angular response to the reaction sphere
	vec3 _magicRotAxis;

	// this axis is that of the rotation matrices applied to the structure
	vec3 _rotationAxis;
	vec3 _rotationCentre;
	vec3 _sphereDiffOffset;
	double _transExponent;
	double _rotExponent;

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
