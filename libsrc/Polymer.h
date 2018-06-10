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
	
	virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal);
	virtual void graph(std::string graphName);
	virtual void differenceGraphs(std::string graphName, CrystalPtr diffCryst);

	static double getBackboneDampening(void *object);
	static void setBackboneDampening(void *object, double value);

	static double getBackboneKick(void *object);
	static void setBackboneKick(void *object, double value);

	static double getSidechainDampening(void *object);
	static void setSidechainDampening(void *object, double value);

	static void setInitialKick(void *object, double value);
	static double getInitialKick(void *object);

	static double getSideKick(void *object);
	static void setSideKick(void *object, double value);

	static double findOverallKickAndDampen(void *object);
	static double vsFindKickAndDampen(void *object);
	
	static double vsSandbox(void *object);
	
	void scaleSidechainsToBFactor();
	void refineBackbone();
	static double vsRefineBackbone(void *object);
	void superimpose();
	
	void attachTargetToRefinement(RefinementStrategyPtr strategy,
	                              FlexGlobal &target);

	static double vsSuperimpose(void *object)
	{
		Parser *parser = static_cast<Parser *>(object);
		Polymer *polymer = dynamic_cast<Polymer *>(parser);
		
		polymer->superimpose();
		return 0;
	}
	
	virtual void reportParameters();
	void downWeightResidues(int start, int end, double value);

	void applyPolymerChanges();
	void refineToEnd(int monNum, CrystalPtr target, RefinementType rType);
	double refineRange(int start, int end, 
	                 CrystalPtr target, RefinementType rType);
	bool test();
	ModelPtr getAnchorModel();
	void findAnchorNearestCentroid();
	void weightStrands();
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

	static void setOverallScale(void *object, double value)
	{
		static_cast<Polymer *>(object)->_overallScale = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}
	
	static double getOverallScale(void *object)
	{
		return static_cast<Polymer *>(object)->_overallScale;
	}
	
	static void vsTransTensorOverall(void *object, double value);
	static double vsFitRotation(void *object);
	static double vsFitTranslation(void *object);

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

	static void vsSetRotAngle(void *object, double value)
	{
		Parser *parser = static_cast<Parser *>(object);
		Polymer *polymer = dynamic_cast<Polymer *>(parser);
		setRotAngle(polymer, value);
		std::cout << "Rot: " << polymer->_rotationAngle << std::endl;
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
	
	static void setRotExponent(void *object, double value)
	{
		Polymer *polymer = static_cast<Polymer *>(object);
		polymer->_rotExponent = value;
		polymer->setChangedRotation();
		polymer->getExtraRotations();
	}

	static double getRotExponent(void *object)
	{
		return static_cast<Polymer *>(object)->_rotExponent;
	}

	static void setTransExponent(void *object, double value)
	{
		Polymer *polymer = static_cast<Polymer *>(object);
		polymer->_transExponent = value;
		polymer->applyTranslationTensor();
	}

	static double getTransExponent(void *object)
	{
		return static_cast<Polymer *>(object)->_transExponent;
	}

	static double getSphereDiffOffsetX(void *object)
	{
		return static_cast<Polymer *>(object)->_sphereDiffOffset.x;
	}

	static void setSphereDiffOffsetX(void *object, double value)
	{
		static_cast<Polymer *>(object)->_sphereDiffOffset.x = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static double getSphereDiffOffsetY(void *object)
	{
		return static_cast<Polymer *>(object)->_sphereDiffOffset.y;
	}

	static void setSphereDiffOffsetY(void *object, double value)
	{
		static_cast<Polymer *>(object)->_sphereDiffOffset.y = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
	}

	static double getSphereDiffOffsetZ(void *object)
	{
		return static_cast<Polymer *>(object)->_sphereDiffOffset.z;
	}

	static void setSphereDiffOffsetZ(void *object, double value)
	{
		static_cast<Polymer *>(object)->_sphereDiffOffset.z = value;
		static_cast<Polymer *>(object)->applyTranslationTensor();
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
	void refineEverything(int start);
	void refineLoop(int start, bool magic);

	mat3x3 _transTensor;
	double _overallScale;
	int _anchorNum;
	double _tmpPhi;
	double _tmpPsi;
	double _startB;
	double _dampening;
	double _kick;
	double _sideDampening;
	double _sideKick;
	int _totalMonomers;
	void minimiseCentroids();
	void minimiseRotations();
	void applyTranslationTensor();

	AtomGroupPtr _allBackbones;

};

#endif /* defined(__vagabond__Polymer__) */
