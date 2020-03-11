//
//  Absolute.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Absolute__
#define __vagabond__Absolute__

#include <stdio.h>
#include "Model.h"
#include "vec3.h"
#include <string>
#include "FileReader.h"
#include "Bond.h"
#include "Parser.h"

/**
 * \class Absolute
 *
 * \brief This class directly adopts the properties describing the
 * flexibility of an atom in a PDB file, and so is a type of Model which
 * explicitly describes an atom in terms of a position (x, y, z) and a B
 * factor (one B or a tensor).
 *
 * These are usually created directly from information from the PDB file.
 * This is best initialised with Absolute(vec3, double, std::string, double).
 *
 * Some Absolute atoms (namely anchors) have had setAnchorPoint() called. This
 * means that this atom can provide a sampling of atom positions, which are
 * namely used as the origin of a recursive bonded structure. The anisotropic
 * tensor is adjusted to compensate for this, and therefore the B factor and
 * anisotropic tensor are not necessarily the same as in the original PDB
 * file.
 */

class Absolute : public Model
{
public:
	/**
	* Best initialiser for Absolute which prepares the basic information it
	* needs.
	* 	\param pos vec3 containing x,y,z of atom positions in Angstroms.
	* 	\param bFac Isotropic B factor describing atom flexibility
	* 	\param element string containing element symbol in upper case.
	* 	\param occValue occupancy between 0 and 1.
	*/
	Absolute(vec3 pos, double bFac, std::string element, double occValue);
	
	/**
	* 	Creates an empty Absolute object and initialises default variables. */
	Absolute();
	virtual ~Absolute() {};

	virtual vec3 getAbsolutePosition()
	{
		return _position;
	}

	/** Makes the corresponding atom when added to the molecule. */
	virtual void addToMonomer(MonomerPtr monomer);
	

	/** Creates an atom from the Absolute model and adds it to the molecule
	 * of interest. Do not use if the Atom already exists. In this case,
	 * add the Atom directly to Molecule using Molecule::addAtom */
	virtual void addToMolecule(MoleculePtr molecule);
	virtual mat3x3 getRealSpaceTensor();
	
	/** Sets a number of fields which would be present in the PDB file. */
	void setIdentity(int resNumValue, std::string chainID,
	                 std::string resName, std::string atomName, int atomNum)
	{
		_resNum = resNumValue;
		_chainID = chainID;
		_resName = resName;
		_atomName = atomName;
		_atomNum = atomNum;

		trim(_chainID); trim(_atomName);

		trim(_resName); to_lower(_resName);
	}

	virtual void propagateChange(int depth, bool refresh)
	{
		Model::propagateChange(depth, refresh);
	}

	void setAlternativeConformerName(std::string conformer)
	{
		_conformer = conformer;
	}

	/** Returns the number of the (amino acid) residue */
	int getResNum()
	{
		return _resNum;
	}
 
	/** Returns the name of the residue, e.g. ser, gly */
	std::string getResName()
	{
		return _resName;
	}

	std::string getChainID()
	{
		return _chainID + (_hetatm ? "_hetatm" : "");
	}

	void setBFactor(double bfac)
	{
		_bFactor = bfac;
		_tensor = make_mat3x3();
		_tensor.vals[0] = bfac;
		_tensor.vals[4] = bfac;
		_tensor.vals[8] = bfac;
	}

	/** Set a tensor to describe an anisotropic B factor. If called, will
	 * return true for isUsingTensor(). */
	void setTensor(mat3x3 tensor);

	double getBFactor()
	{
		return _bFactor;
	}

	void setHeteroAtom(bool hetatm)
	{
		_hetatm = hetatm;
	}

	bool isHeteroAtom()
	{
		return _hetatm;
	}
	
	/** Returns if the Absolute model uses an anisotropic B factor for whatever
	 * reason. Certainly the case if loaded from a PDB file with an ANISO
	 * record. */
	bool isUsingTensor()
	{
		return _usingTensor;	
	}

	virtual double getMeanSquareDeviation();

	static double getB(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_bFactor;
	}

	static void setB(void *object, double b)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_bFactor = b;
	}

	virtual std::string getClassName()
	{
		return "Absolute";
	}

	/** Returns the atom which this model directly controls. */
	virtual AtomPtr getAtom()
	{
		return _atom.lock();
	}
	
	/** Nothing to be done in this class. */
	void refreshPositions() {};
	
	virtual bool hasExplicitPositions()
	{
		return false;
	}

	void resetSamples();
	
	virtual double getEffectiveOccupancy()
	{
		return _occupancy;	
	}
	
	vec3 getRandomPosition();

	static void setPosX(void *object, double x)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_position.x = x;
	}

	static void setPosY(void *object, double y)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_position.y = y;
	}

	static void setPosZ(void *object, double z)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		abs->_position.z = z;
	}

	static double getPosX(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_position.x;
	}

	static double getPosY(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_position.y;
	}

	static double getPosZ(void *object)
	{
		Absolute *abs = static_cast<Absolute *>(object);
		return abs->_position.z;
	}
	
	std::string getElementSymbol()
	{
		return _element;
	}
	
	AtomPtr makeAtom();
	
	virtual std::string shortDesc();
protected:
	static double getExpValue(void *object, double x, double y, double z);

	virtual std::string getParserIdentifier()
	{
		return "absolute_" + i_to_str(getAtom()->getAtomNum());
	}

	virtual void addProperties();
	virtual void addObject(ParserPtr, std::string) {};
	virtual void linkReference(BaseParserPtr object, std::string category);
	virtual void postParseTidy();
private:
	void initialise();
	AtomWkr _atom;
	std::vector<AtomWkr> _nextAtoms;
	std::string _element;
	double _occupancy;
	std::string _chainID, _resName, _atomName;
	int _resNum, _atomNum;
	bool _hetatm;
	bool _usingTensor;
	mat3x3 _tensor;
	std::string _conformer;
	std::vector<BondSample> _bondSamples;

	vec3 _position;
	double _bFactor;
};

#endif /* defined(__vagabond__Absolute__) */
