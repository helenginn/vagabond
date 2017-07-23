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

class Atom : public std::enable_shared_from_this<Atom>
{
public:
	Atom();

	void setModel(ModelPtr model);

	bool isBackbone();
	bool isBackboneAndSidechain();

	vec3 getPosition()
	{
		return _position;
	}

	void setPosition(vec3 pos)
	{
		_position = pos;
	}

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

	/* Returns a FFT for the model dist, for reuse */
	void addToMap(FFTPtr fft, mat3x3 unit_cell);

private:
	ModelPtr _model;
	ElementPtr _element;
	std::string _atomName;
	MonomerWkr _monomer;

	vec3 _position;
};

#endif /* defined(__vagabond__Atom__) */
