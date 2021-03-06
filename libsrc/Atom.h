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
#include "Element.h"
#include <string>
#include <mutex>
#include "../libinfo/GeomTable.h"
#include "Parser.h"
#include "FileReader.h"

/**
 * \class Atom
 * \brief Responsible for the PDB-derived properties of a
 * single ATOM or HETATM line from a PDB.
 */

class Atom : public Parser
{
public:
	AtomPtr shared_from_this()
	{
		return ToAtomPtr(Parser::shared_from_this());
	}
	
	Atom();
	Atom(Atom &other);

	virtual ~Atom();

	/** Change the atom's model to a new one, which will now be called when
	 * the atom distribution is required. Can be changed at any time */
	void setModel(ModelPtr model);

	/** Atom is *only* part of the backbone (so *not* including C-alpha atoms)
	 */
	bool isBackbone();
	
	/** If atom is considered both part of backbone and sidechain, such as
	* C-alphas for protein chains. */
	bool isBackboneAndSidechain();

	bool isSidechain()
	{
		return !(isBackbone() || isBackboneAndSidechain());
	}

	/** Get the average absolute position from the atom's Model */
	virtual vec3 getAbsolutePosition();
	
	/** Positional displacement between average ensemble position and
	* reference position (usually from PDB). */
	double posDisplacement(bool fromSaved = false, bool refresh = true,
	                       bool sq = false);
	
	void saveInitialPosition();

	vec3 getPositionInUnitCell();

	/** Change the periodic table element of the atom */
	void setElement(ElementPtr element)
	{
		_element = element;
		_elementSymbol = element->getSymbol();
	}
	
	/** Set the monomer for this atom with no other connections. Do not use if
	* you've made a new atom to assign to a monomer. Use Monomer::addAtom()
	* instead. */
	void setMonomer(MonomerPtr monomer)
	{
		_monomer = monomer;
	}

	/** Change atom name such as CA, CB etc. */
	void setAtomName(std::string name)
	{
		_atomName = name;
	}

	/** Get the atom name, such as CA, CB etc. */
	std::string getAtomName()
	{
		return _atomName;
	}

	/** Get the model associated with the position/flexibility of this atom */
	ModelPtr getModel()
	{
		return _model;
	}
	
	/** Get the model associated with the position/flexibility of this atom
	 * 	and cast as an explicit model (only if you're certan it is! */
	ExplicitModelPtr getExplicitModel();
	
	/** Call the appropriate element's Element::electronCount() - convenience
	 * function */
	int getElectronCount();

	/** Get the associated periodic table element for this Atom */
	ElementPtr getElement()
	{
		return _element;
	}

	/** Get the radius of the atom used for solvent calculations */
	double getSolventRadius();

	/** Solvent content added for multiple conformers of a mask at once,
	 * with one flag per bit */
	void addManyToMask(VagFFTPtr fft, int conf, int total);
	void addToSolventMask(VagFFTPtr fft, double radius,
	                      std::vector<Atom *> *ptrs, int conf);
	void addPointerToLocalArea(VagFFTPtr fft, vec3 pos,
	                           std::vector<Atom *> *ptrs,
	                           double rad = 0);

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

	virtual vec3 getInitialPosition()
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
	void writePositionsToFile(std::string suffix = "");

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
	virtual std::string shortDesc();
	virtual std::string longDesc();

	MoleculePtr getMolecule();
	std::string getPDBContribution(int ensembleNum = -1);
	std::string averagePDBContribution(bool samePos, bool sameB);
	std::string anisouPDBLine(CrystalPtr crystal);

	/* Tolerance in Angstroms. */
	bool closeToAtom(AtomPtr another, double tolerance = 2);
	
	void cacheCloseWaters(double tolerance = 5);

	double getDistanceFrom(Atom *other, int nSample = -1, bool quick = false);
	static double getAngle(AtomPtr atom1, AtomPtr atom2, AtomPtr atom3);
	
	std::string getChainID();
	
	void setWater(int set = 1)
	{
		_isWater = set;
	}

	double posToMouse();
	
	void setTargetPosition(vec3 pos, double weight)
	{
		_targetPos = pos;
		_targetWeight = weight;
	}
	
	bool isAtom(std::string atomName, int resNum);
	
	double getTargetWeight()
	{
		return _targetWeight;
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
	
	void clearTargetB()
	{
		_targetB = 0;
		_targetBCount = 0;
	}
	
	double getTargetB()
	{
		return _targetB / _targetBCount;
	}
	
	void setTargetB(double b)
	{
		_targetB += b;
		_targetBCount++;
	}
	
	void setWeightOnly(double mult)
	{
		_weightOnly = mult;
	}
	
	vec3 getPositionInAsu(int conf = -1);
	
	void setGhostBond(GhostBondPtr ghost)
	{
		_ghost = ghost;
	}
	
	GhostBondPtr getGhostBond()
	{
		return _ghost;
	}
	
	std::string getElementSymbol()
	{
		return _elementSymbol;
	}

	double fishWhackMagnitude();
	
	static bool greater(AtomPtr a1, AtomPtr a2);
	vec3 getSymRelatedPosition(int i, vec3 pos);
protected:
	vec3 getSymRelatedPosition(int i, int conf = -1);

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
	size_t symOpCount();
	ModelPtr _model;
	ModelPtr _distModelOnly;
	ElementPtr _element;
	std::string _atomName;
	MonomerWkr _monomer;
	vec3 _initialPosition;
	double _initialB;
	double _targetB;
	int _targetBCount;
	vec3 _pdbPosition;
	int _atomNum;
	int _asu;
	double _origOccupancy;
	vec3 _ellipsoidLongestAxis;
	double _weighting;
	double _weightOnly;
	std::string _conformer;
	std::string _elementSymbol;
	bool _fromPDB;
	int _isWater;
	int _hetatm;
	mat3x3 _tensor;
	bool _hBondage;

	vec3 _targetPos;
	double _targetWeight;

	GhostBondPtr _ghost;
	AtomType _geomType;

	std::mutex _whackLock;
};

#endif /* defined(__vagabond__Atom__) */
