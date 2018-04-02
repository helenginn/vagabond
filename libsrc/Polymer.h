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

class Polymer : public Molecule
{
public:
	Polymer()
	{
		_dampening = 0.05;
		_sideDampening = 0.05;
		_sideKick = 0;
		_anchorNum = 0;
		_totalMonomers = 0;
		_transTensor = make_mat3x3();
		_startB = Options::getBStart();
		_extraRotParams = {1, 0, 0};
		_tmpPhi = 0;
		_tmpPsi = 0;
	}

	void closenessSummary();
	void addMonomer(MonomerPtr monomer);
	virtual void summary();
	virtual void tieAtomsUp();
	void splitConformers();
	virtual void refine(CrystalPtr target, RefinementType rType);
	virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal);
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

	void scaleSidechainsToBFactor();
	void superimpose();
	virtual void reportParameters();
	void downWeightResidues(int start, int end, double value);

	void applyPolymerChanges();
	void refineToEnd(int monNum, CrystalPtr target, RefinementType rType);
	void refineRange(int start, int end, 
	                 CrystalPtr target, RefinementType rType);
	bool test();
	ModelPtr getAnchorModel();
	void findAnchorNearestCentroid();
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

	long monomerCount()
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

	static double getTransTensor11(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[0];
	}

	static double getTransTensor21(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[3];
	}

	static double getTransTensor12(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[1];
	}

	static double getTransTensor31(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[6];
	}

	static double getTransTensor13(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[2];
	}

	static double getTransTensor22(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[4];
	}

	static double getTransTensor23(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[5];
	}

	static double getTransTensor32(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[7];
	}

	static double getTransTensor33(void *object)
	{
		return static_cast<Polymer *>(object)->_transTensor.vals[8];
	}

	static void setTransTensor11(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[0] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor12(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[1] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor21(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[3] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor13(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[2] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor31(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[6] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor22(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[4] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor32(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[7] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor23(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[5] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setTransTensor33(void *object, double value)
	{
		static_cast<Polymer *>(object)->_transTensor.vals[8] = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static void setRotPhi(void *object, double value)
	{
		Polymer *polymer = static_cast<Polymer *>(object);
		polymer->_tmpPhi = value;
		polymer->refreshRotationAxis();
		polymer->getExtraRotations();
	}

	static double getRotPhi(void *object)
	{
		return static_cast<Polymer *>(object)->_tmpPhi;
	}

	static void setRotPsi(void *object, double value)
	{
		Polymer *polymer = static_cast<Polymer *>(object);
		polymer->_tmpPsi = value;
		polymer->refreshRotationAxis();
		polymer->getExtraRotations();
	}

	static double getRotPsi(void *object)
	{
		return static_cast<Polymer *>(object)->_tmpPsi;
	}

	static void setRotAngle(void *object, double value)
	{
		Polymer *polymer = static_cast<Polymer *>(object);
		polymer->_rotationAngle = value;
		polymer->setChangedRotation();
		polymer->getExtraRotations();
	}

	static double getRotAngle(void *object)
	{
		return static_cast<Polymer *>(object)->_rotationAngle;
	}

	virtual void calculateExtraRotations();
	std::vector<vec3> getAnchorSphereDiffs();
	void optimiseWholeMolecule(bool translation = false, bool rotation = true);
	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void postParseTidy();
	
	AtomGroupPtr getAllBackbone();
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
	vec3 _extraRotParams;

	void refreshRotationAxis()
	{
		mat3x3 rot = mat3x3_rot_from_angles(_tmpPhi, _tmpPsi);
		_rotationAxis = mat3x3_axis(rot, 0);
		setChangedRotation();
	}

	mat3x3 _transTensor;
	int _anchorNum;
	double _tmpPhi;
	double _tmpPsi;
	double _startB;
	double _dampening;
	double _sideDampening;
	double _sideKick;
	double _totalMonomers;
	void minimiseCentroids();
	void minimiseRotations();
	void applyTranslationTensor();

	AtomGroupPtr _allBackbones;

};

#endif /* defined(__vagabond__Polymer__) */
