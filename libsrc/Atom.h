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
	vec3 getAbsolutePosition();
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
	std::string pdbLineBeginning(std::string start = "ATOM  ");

	AtomType getGeomType()
	{
		return _geomType;
	}

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

	void setEllipsoidLongestAxis(vec3 axis)
	{
		_ellipsoidLongestAxis = axis;
	}

	vec3 getEllipsoidLongestAxis()
	{
		return _ellipsoidLongestAxis;
	}

	std::string shortDesc();

	MoleculePtr getMolecule();
	void setKeepModel();
	std::string getPDBContribution(int ensembleNum = -1);
	std::string averagePDBContribution(bool samePos, bool sameB);
	std::string anisouPDBLine(CrystalPtr crystal);

	static double getAngle(AtomPtr atom1, AtomPtr atom2, AtomPtr atom3);
	
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
	vec3 _ellipsoidLongestAxis;
	double _weighting;

	AtomType _geomType;
};

#endif /* defined(__vagabond__Atom__) */
