//
//  Atom.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Atom__
#define __vagabond__Atom__

#include <stdio.h>
#include "shared_ptrs.h"
#include <vector>
#include "vec3.h"
#include "mat3x3.h"
#include <string>
#include "../libinfo/GeomTable.h"
#include "fftw3d.h"
#include "Parser.h"
#include "FileReader.h"

/**
 * \class Atom
 * \brief Responsible for the PDB-derived properties of a
 * single ATOM or HETATM line from a PDB.
 */

class Plucker;

class Atom : public Parser
{
public:
	AtomPtr shared_from_this()
	{
		return ToAtomPtr(Parser::shared_from_this());
	}
	
	Atom();
	Atom(Atom &other);

	virtual ~Atom() {}	

	void setModel(ModelPtr model);
	FFTPtr getBlur();

	bool isBackbone();
	
	/** If atom is considered both part of backbone and sidechain, such as
	* C-alphas for protein chains. */
	bool isBackboneAndSidechain();

	vec3 getAbsolutePosition();
	vec3 getAsymUnitPosition(CrystalPtr crystal, int nSample = -1);
	
	/** Positional displacement between average ensemble position and
	* reference position (usually from PDB). */
	double posDisplacement();

	void setElement(ElementPtr element)
	{
		_element = element;
	}
	
	/** Set the monomer for this atom with no frilly bits. Do not use if
	* you've made a new atom to assign to a monomer. Use Monomer::addAtom()
	* instead. */
	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}

	void setAtomName(std::string name)
	{
		_atomName = name;
	}

	std::string getAtomName()
	{
		return _atomName;
	}

	ModelPtr getModel()
	{
		return _model;
	}
	
	ExplicitModelPtr getExplicitModel();
	
	int getElectronCount();

	ElementPtr getElement()
	{
		return _element;
	}

	double getSolventRadius();
	void addToSolventMask(FFTPtr fft, mat3x3 unit_cell, double radius = -1);
	void addPointerToLocalArea(FFTPtr fft, mat3x3 unit_cell,
	                           std::vector<Atom *> *ptrs,
	                           double rad = 0);

	/* Returns a FFT for the model dist, for reuse */
	void addToMap(FFTPtr fft, mat3x3 unit_cell,
	              vec3 offset = make_vec3(0, 0, 0), bool mask = false,
	              bool sameScale = false, bool noWrap = false);

	void setOriginalOccupancy(double occ)
	{
		_origOccupancy = occ;
	}

	double getOriginalOccupancy()
	{
		return _origOccupancy;
	}

	void setAlternativeConformer(std::string conf)
	{
		_conformer = conf;
	}

	std::string getAlternativeConformer()
	{
		return _conformer;
	}

	bool isFromPDB()
	{
		return _fromPDB;
	}

	void setFromPDB(bool value)
	{
		_fromPDB = value;
	}

	void setInitialPosition(vec3 pos)
	{

		_initialPosition = pos;
	}

	vec3 getInitialPosition()
	{
		return _initialPosition;
	}

	vec3 getPDBPosition()
	{
		return _pdbPosition;
	}

	void setPDBPosition(vec3 pdbPos)
	{
		_pdbPosition = pdbPos;
	}
	
	double getBFactor();

	double getInitialBFactor()
	{
		return _initialB;
	}

	void setInitialBFactor(double b)
	{
		_initialB = b;
		double inAng = b / (8 * M_PI * M_PI);
		_tensor = make_mat3x3();
		mat3x3_mult_scalar(&_tensor, inAng);
	}

	void setAtomNum(int atomNum)
	{
		_atomNum = atomNum;
	}

	int getAtomNum()
	{
		return _atomNum;
	}

	void findAtomType(std::string resName);
	void inheritParents();
	std::string pdbLineBeginning(std::string start = "ATOM  ");
	void writePositionsToFile();

	AtomType getGeomType();
	void convertToDisulphide();
	void refreshBondAngles();

	MonomerPtr getMonomer()
	{
		if (_monomer.expired())
		{
			return MonomerPtr();
		}

		return _monomer.lock();
	}

	double getWeighting()
	{
		return _weighting;
	}

	void setWeighting(double weighting)
	{
		_weighting = weighting;
	}

	void setHetatm(int hetatm)
	{
		_hetatm = hetatm;
	}
	
	bool isHeteroAtom()
	{
		return (_hetatm == 1);
	}

	void setEllipsoidLongestAxis(vec3 axis)
	{
		_ellipsoidLongestAxis = axis;
	}

	vec3 getEllipsoidLongestAxis()
	{
		return _ellipsoidLongestAxis;
	}

	void setTensor(mat3x3 tensor)
	{
		_tensor = tensor;
	}

	/**
	* Get the tensor as described in the original PDB file.
	* This may be different from the tensor rederived from Vagabond models.
	* \return symmetric matrix with 6 unique values
	*/
	mat3x3 getTensor()
	{
		return _tensor;
	}
	
	int getResidueNum();

	std::string description();
	std::string shortDesc();

	MoleculePtr getMolecule();
	std::string getPDBContribution(int ensembleNum = -1);
	std::string averagePDBContribution(bool samePos, bool sameB);
	std::string anisouPDBLine(CrystalPtr crystal);

	/* Tolerance in Angstroms. */
	bool closeToAtom(AtomPtr another, double tolerance = 2);
	
	void cacheCloseWaters(double tolerance = 5);

	double getDistanceFrom(Atom *other, int nSample = -1, bool quick = false);
	static double getAngle(AtomPtr atom1, AtomPtr atom2, AtomPtr atom3);
	
	size_t pluckCount();
	
	void setWater(int set = 1)
	{
		_isWater = set;
	}
	
	bool isWater()
	{
		return _isWater;
	}
	
	bool canBeHydrogenBonder()
	{
		return (_hBondage);
	}
	
	void setHBonding(bool status)
	{
		_hBondage = status;
	}
	
	Atom *pluckAnother();
protected:
	virtual std::string getClassName()
	{
		return "Atom";
	}

	virtual std::string getParserIdentifier()
	{
		return "atom_" + i_to_str(_atomNum) + "_" + shortDesc(); 
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void postParseTidy();
private:
	vec3 getPositionInAsu();
	vec3 getSymRelatedPosition(int i);
	size_t symOpCount();
	ModelPtr _model;
	ModelPtr _distModelOnly;
	ElementPtr _element;
	std::string _atomName;
	MonomerWkr _monomer;
	vec3 _initialPosition;
	double _initialB;
	vec3 _pdbPosition;
	int _atomNum;
	int _asu;
	double _origOccupancy;
	vec3 _ellipsoidLongestAxis;
	double _weighting;
	std::string _conformer;
	std::string _elementSymbol;
	bool _fromPDB;
	int _isWater;
	int _hetatm;
	mat3x3 _tensor;
	bool _hBondage;

	AtomType _geomType;
	
	Plucker *_waterPlucker;
};

#endif /* defined(__vagabond__Atom__) */
