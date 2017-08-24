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

class Atom : public std::enable_shared_from_this<Atom>
{
public:
	Atom();
	Atom(Atom &other);

	void setModel(ModelPtr model);
	FFTPtr getBlur();

	bool isBackbone();
	bool isBackboneAndSidechain();

	vec3 getPosition();

	void setElement(ElementPtr element)
	{
		_element = element;
	}

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

	ElementPtr getElement()
	{
		return _element;
	}

	/* Fit with FFT for the element dist only */
	double scoreWithMap(FFTPtr fft, mat3x3 unit_cell,
						std::vector<double> *xs = NULL,
						std::vector<double> *ys = NULL);

	/* Returns a FFT for the model dist, for reuse */
	void addToMap(FFTPtr fft, mat3x3 unit_cell, vec3 offset = make_vec3(0, 0, 0));

	void setInitialPosition(vec3 pos)
	{
		_initialPosition = pos;
	}

	vec3 getInitialPosition()
	{
		return _initialPosition;
	}

	double getInitialBFactor()
	{
		return _initialB;
	}

	void setInitialBFactor(double b)
	{
		_initialB = b;
	}

	void setAtomNum(int atomNum)
	{
		_atomNum = atomNum;
	}

	void findAtomType(std::string resName);
	void inheritParents();
	std::string pdbLineBeginning(int i);

	AtomType getGeomType()
	{
		return _geomType;
	}

	MonomerPtr getMonomer()
	{
		return _monomer.lock();
	}

	void setKeepModel();
private:
	ModelPtr _model;
	ModelPtr _distModelOnly;
	ElementPtr _element;
	std::string _atomName;
	MonomerWkr _monomer;
	vec3 _initialPosition;
	double _initialB;
	int _atomNum;

	AtomType _geomType;
};

#endif /* defined(__vagabond__Atom__) */
