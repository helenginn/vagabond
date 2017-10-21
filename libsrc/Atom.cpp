//
//  Atom.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Atom.h"
#include "fftw3d.h"
#include "mat3x3.h"
#include <math.h>
#include <stdlib.h>
#include "Bond.h"
#include <iostream>
#include "Element.h"
#include "Monomer.h"
#include "Polymer.h"
#include <iomanip>
#include "FileReader.h"
#include <sstream>
#include "Crystal.h"
#include "PDBReader.h"

Atom::Atom()
{
	_initialPosition = make_vec3(0, 0, 0);
	_initialB = 0;
	_geomType = AtomUnassigned;
}

Atom::Atom(Atom &other)
{
	_initialPosition = other._initialPosition;
	_initialB = other._initialB;
	_geomType = other._geomType;
	_element = other._element;
	_atomName = other._atomName;
	_model = other._model;
	_monomer = other._monomer;
}

void Atom::inheritParents()
{
	getMonomer()->addAtom(shared_from_this());
	getMonomer()->getPolymer()->addAtom(shared_from_this());
}

void Atom::setModel(ModelPtr model)
{
	_model = model;
}

FFTPtr Atom::getBlur()
{
	FFTPtr modelDist = _model->getDistribution();
	/* YEAH, that line does something */

	if (_distModelOnly)
	{
		modelDist = _distModelOnly->getDistribution();
	}

	return modelDist;
}

double Atom::scoreWithMap(CrystalPtr crystal, std::vector<double> *xs,
						  std::vector<double> *ys, bool diff,
						  MapScoreType mapScore)
{
	FFTPtr fft = crystal->getFFT();

	if (diff)
	{
		fft = crystal->getDiFFT();
	}

	return scoreWithMap(fft, crystal->getReal2Frac(), xs, ys, mapScore);
}

double Atom::scoreWithMap(FFTPtr fft, mat3x3 unit_cell,
						  std::vector<double> *xs, std::vector<double> *ys,
						  MapScoreType mapScore)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modelDist = getBlur();
	FFT::multiply(modelDist, atomDist);
	modelDist->fft(1);
	modelDist->invertScale();

	vec3 pos = _model->getAbsolutePosition();
	mat3x3_mult_vec(unit_cell, &pos);

	if (pos.x != pos.x)
	{
		return 0;
	}

	double score = FFT::score(fft, modelDist, pos, xs, ys, mapScore);

	return score;
}

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell, vec3 offset, bool useNew)
{
	FFTPtr atomDist = _element->getDistribution();
	FFTPtr modified;

	if (!useNew)
	{
		modified = getBlur();
	}
	else
	{
		modified = _model->getDistribution();
	}

	FFT::multiply(modified, atomDist);
	modified->fft(1);
	modified->invertScale();

	vec3 pos = _model->getAbsolutePosition();

	if (_distModelOnly)
	{
		pos = _distModelOnly->getAbsolutePosition();
	}

	pos = vec3_subtract_vec3(pos, offset);
	mat3x3_mult_vec(unit_cell, &pos);

	if (pos.x != pos.x)
	{
		return;
	}

	FFT::add(fft, modified, pos);
}

vec3 Atom::getPosition()
{
	return _model->getStaticPosition();
}

bool Atom::isBackbone()
{
	if (_atomName == "C") return true;
	if (_atomName == "H") return true;
	if (_atomName == "N") return true;
	if (_atomName == "O") return true;
	if (_atomName == "HA3") return true;

	return false;
}

bool Atom::isBackboneAndSidechain()
{
	if (_atomName == "CA") return true;
	if (_atomName == "HA") return true;
	if (_atomName == "HA2") return true;

	return false;
}

/* Convert to lookup table */
void Atom::findAtomType(std::string resName)
{
	_geomType = GeomTable::getGeomTable().getType(resName, _atomName);
}

std::string Atom::pdbLineBeginning(int i)
{
	std::string residueName = getMonomer()->getIdentifier();
	int resNum = getMonomer()->getResidueNum();
	to_upper(residueName);
	std::ostringstream line;

	char conformer[] = "A";
	conformer[0] += i;

	line << "ATOM  ";
	line << std::setfill(' ') << std::setw(5) << std::fixed << _atomNum;
	line << std::setfill(' ') << std::setw(4) << _atomName;
	line << " " << conformer;
	line << std::setw(3) << residueName;
	line << " A";
	line << std::setfill(' ') << std::setw(4) << resNum;
	line << "    ";

	return line.str();
}

void Atom::setKeepModel()
{
	_distModelOnly = _model;
}

double Atom::posDisplacement()
{
	BondPtr bond = ToBondPtr(getModel());
	bond->getDistribution();
	vec3 bestPos = bond->getAbsolutePosition();
	vec3 initialPos = getPDBPosition();

	vec3 diff = vec3_subtract_vec3(bestPos, initialPos);
	double score = vec3_length(diff);

	return score;
}

std::string Atom::averagePDBContribution(bool samePos, bool sameB)
{
	getModel()->getDistribution();
	std::string atomName = getAtomName();
	ElementPtr element = getElement();

	double occupancy = 1;
	vec3 placement = getModel()->getAbsolutePosition();

	if (samePos)
	{
		placement = getPDBPosition();
	}

	double bFactor = getModel()->getMeanSquareDeviation();

	if (sameB)
	{
		bFactor = getInitialBFactor();
	}

	return PDBReader::writeLine(shared_from_this(),
								placement, 0, occupancy, bFactor);
}

std::string Atom::getPDBContribution()
{
	std::string atomName = getAtomName();
	ElementPtr element = getElement();

	int tries = 10;

	if (element->getSymbol() == "H")
	{
		tries = 1;
	}

	std::ostringstream stream;
	std::vector<BondSample> positions = getModel()->getFinalPositions();

	double skip = (double)positions.size() / 25.;

	if (skip < 0) skip = 1;
	const int side = 7;
	int count = 0;

	for (double i = 0; i < positions.size(); i+= 1)
	{
		int l = i / (side * side);
		int k = (i - (l * side * side)) / side;
		int h = (i - l * side * side - k * side);

		if ((h + k) % 2 != 0 || (k + l) % 2 != 0 || (l + h) % 2 != 0)
		{
			continue;
		}

		vec3 placement = positions[i].start;
		double occupancy = 50 * positions[i].occupancy / double(tries);
		stream << PDBReader::writeLine(shared_from_this(), placement, count, occupancy, 0);
		count++;
	}

	return stream.str();
}

std::string Atom::shortDesc()
{
	return getMonomer()->getIdentifier()
		+ i_to_str(getMonomer()->getResidueNum()) +
		getAtomName();
}

MoleculePtr Atom::getMolecule()
{
	if (getMonomer())
	{
		return getMonomer()->getPolymer();
	}

	return MoleculePtr();
}
