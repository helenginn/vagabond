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
	double posDisplacement();

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
						std::vector<double> *ys = NULL,
						MapScoreType mapScore = MapScoreTypeCorrel);

	double scoreWithMap(CrystalPtr crystal, std::vector<double> *xs = NULL,
						std::vector<double> *ys = NULL, bool diff = false,
						MapScoreType mapScore = MapScoreTypeCorrel);

	/* Returns a FFT for the model dist, for reuse */
	void addToMap(FFTPtr fft, mat3x3 unit_cell,
				  vec3 offset = make_vec3(0, 0, 0), bool useNew = false);

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

	double getInitialBFactor()
	{
		return _initialB;
	}

	double getInitialAnisoB(int index)
	{
		return _aniso[index];
	}

	void setInitialAnisoBs(double x, double y, double z)
	{
		const double bFacMult = 8 * M_PI * M_PI / 3;
		_aniso[0] = x * bFacMult;
		_aniso[1] = y * bFacMult;
		_aniso[2] = z * bFacMult;
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

	std::string shortDesc();


	void setKeepModel();
	std::string getPDBContribution();
private:
	ModelPtr _model;
	ModelPtr _distModelOnly;
	ElementPtr _element;
	std::string _atomName;
	MonomerWkr _monomer;
	vec3 _initialPosition;
	double _initialB;
	vec3 _pdbPosition;
	int _atomNum;
	double _aniso[3];

	AtomType _geomType;
};

#endif /* defined(__vagabond__Atom__) */
