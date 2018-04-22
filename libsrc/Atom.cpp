//
//  Atom.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Atom.h"
#include "Absolute.h"
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
#include "Anisotropicator.h"

Atom::Atom()
{
	_initialPosition = make_vec3(0, 0, 0);
	_pdbPosition = make_vec3(0, 0, 0);
	_ellipsoidLongestAxis = make_vec3(0, 0, 0);
	_initialB = 0;
	_geomType = AtomUnassigned;
	_weighting = 1;
	_origOccupancy = 1.0;
	_fromPDB = true;
	_tensor = make_mat3x3();
	_hetatm = -1;
}

Atom::Atom(Atom &other)
{
	_initialPosition = other._initialPosition;
	_initialB = other._initialB;
	_geomType = other._geomType;
	_element = other._element;
	_atomName = other._atomName;
	_model = other._model;
	_distModelOnly = other._distModelOnly;
	_monomer = other._monomer;
	_weighting = other._weighting;
	_atomNum = other._atomNum;
	_fromPDB = other._fromPDB;
	_pdbPosition = other._pdbPosition;
	_origOccupancy = other._origOccupancy;
	_conformer = other._conformer;
	_ellipsoidLongestAxis = other._ellipsoidLongestAxis;
}

int Atom::getResidueNum()
{
	if (getMonomer())
	{
		return getMonomer()->getResidueNum();	
	}	
	
	return 0;
}

AtomType Atom::getGeomType()
{
	if (_geomType == AtomUnassigned)
	{
		std::string id = getMonomer()->getIdentifier();
		findAtomType(id);
	}

	return _geomType;
}

void Atom::convertToDisulphide()
{
	if (getGeomType() == AtomCysCB)
	{
		_geomType = AtomCysCBS;
	}	
	else if (getGeomType() == AtomCysSG)
	{
		_geomType = AtomCysSGS;	
	}
	else
	{
		return;
	}
	
	if (!getModel()->isBond())
	{
		return;
	}
	
	BondPtr bond = ToBondPtr(getModel());
	ModelPtr parent = bond->getParentModel();
	
	if (!parent->isBond())
	{
		return;	
	}
	
	ToBondPtr(parent)->resetBondAngles();
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

	modelDist->multiplyAll(_weighting);

	return modelDist;
}

double Atom::scoreWithMap(CrystalPtr crystal, std::vector<CoordVal> *vals,
                          bool diff, MapScoreType mapScore)
{
	FFTPtr fft = crystal->getFFT();

	if (diff)
	{
		fft = crystal->getDiFFT();
	}

	return scoreWithMap(fft, crystal->getReal2Frac(), vals, mapScore);
}

double Atom::scoreWithMap(FFTPtr fft, mat3x3 unit_cell,
                          std::vector<CoordVal> *vals,
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

	double score = FFT::score(fft, modelDist, pos, vals, mapScore);

	return score;
}

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell, vec3 offset, bool mask)
{
	FFTPtr atomDist, modified;
	
	if (getModel()->getEffectiveOccupancy() <= 0) return;
	
	if (!mask)
	{
		modified = getBlur();
		atomDist = _element->getDistribution(false, modified->nx);
		FFT::multiply(modified, atomDist);
		modified->fft(1);
		modified->invertScale();
		modified->setTotal(_element->electronCount() * 10e4);
		double occ = _model->getEffectiveOccupancy();
		modified->multiplyAll(occ);
	}
	else
	{
		modified = _element->getMask();		
		double occ = _model->getEffectiveOccupancy();
		modified->multiplyAll(occ);
	}

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

vec3 Atom::getAbsolutePosition()
{
	return _model->getAbsolutePosition();
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

std::string Atom::pdbLineBeginning(std::string start)
{
	std::string residueName = getMonomer()->getIdentifier();
	int resNum = getMonomer()->getResidueNum();
	to_upper(residueName);
	std::ostringstream line;

	char conformer[] = " ";
	if (getAlternativeConformer().length())
	{
		conformer[0] = getAlternativeConformer()[0];
	}

	line << start;
	line << std::setfill(' ') << std::setw(5) << std::fixed << _atomNum;
	line << " " << std::right << std::setfill(' ') << std::setw(4) << _atomName;
	line << std::right << conformer;
	line << std::setw(3) << residueName;
	line << " " << getMolecule()->getChainID()[0];
	line << std::setfill(' ') << std::setw(4) << resNum;
	line << "  ";

	return line.str();
}

double Atom::fullPositionDisplacement()
{
	if (!isFromPDB())
	{
		return 0;
	}

	std::vector<BondSample> samples = getModel()->getFinalPositions();
	vec3 initialPos = getInitialPosition();
	initialPos = getPDBPosition();

	double score = 0;
	
	for (size_t i = 0; i < samples.size(); i++)
	{
		vec3 aPos = samples[i].start;
		vec3 diff = vec3_subtract_vec3(aPos, initialPos);
		double sqlength = vec3_sqlength(diff);
		score += sqlength;
	}
	
	score /= (double)samples.size();

	return score;
}

double Atom::posDisplacement()
{
	if (!isFromPDB())
	{
		return 0;
	}

	getModel()->getFinalPositions();
	vec3 bestPos = getModel()->getAbsolutePosition();
	vec3 initialPos = getPDBPosition();

	vec3 diff = vec3_subtract_vec3(bestPos, initialPos);
	double score = vec3_length(diff);

	return score;
}

std::string Atom::anisouPDBLine(CrystalPtr)
{
	if (!getMonomer() || getAlternativeConformer().length())
	{
		return "";
	}

	std::ostringstream stream;
	stream << pdbLineBeginning("ANISOU");

	mat3x3 realTensor = getModel()->getRealSpaceTensor();
	mat3x3 recipTensor = realTensor;

	double scale = 10e4;
	stream << std::setprecision(0) << std::fixed;
	stream << std::setfill(' ') << std::setw(7) << scale * recipTensor.vals[0];
	stream << std::setfill(' ') << std::setw(7) << scale * recipTensor.vals[4];
	stream << std::setfill(' ') << std::setw(7) << scale * recipTensor.vals[8];
	stream << std::setfill(' ') << std::setw(7) << scale * recipTensor.vals[1];
	stream << std::setfill(' ') << std::setw(7) << scale * recipTensor.vals[2];
	stream << std::setfill(' ') << std::setw(7) << scale * recipTensor.vals[5];
	stream << std::endl;

	return stream.str();
}

std::string Atom::averagePDBContribution(bool samePos, bool sameB)
{
	if (!getMonomer())
	{
		return "";
	}

	getModel()->getFinalPositions();
	std::string atomName = getAtomName();
	ElementPtr element = getElement();

	double occupancy = 1;
	occupancy = getModel()->getEffectiveOccupancy();

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

std::string Atom::getPDBContribution(int ensembleNum)
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

		if ((ensembleNum < 0) &&
		    ((h + k) % 2 != 0 || (k + l) % 2 != 0 || (l + h) % 2 != 0))
		{
			continue;
		}

		if (ensembleNum >= 0)
		{
			i = ensembleNum;
		}

		vec3 placement = positions[i].start;
		double occupancy = 50 * positions[i].occupancy / double(tries);
		stream << PDBReader::writeLine(shared_from_this(), placement, count, occupancy, 0);

		if (ensembleNum >= 0)
		{
			break;
		}

		count++;
	}

	return stream.str();
}

std::string Atom::shortDesc()
{
	if (!getMonomer())
	{
		return getAtomName() + "_" + _conformer;
	}

	std::string start = getMonomer()->getIdentifier()
	+ i_to_str(getMonomer()->getResidueNum())
	+ getAtomName();

	if (_conformer.length())
	{
		start += "_" + _conformer;
	}

	return start;
}

MoleculePtr Atom::getMolecule()
{
	if (getMonomer())
	{
		return getMonomer()->getPolymer();
	}

	return MoleculePtr();
}

double Atom::getAngle(AtomPtr atom1, AtomPtr atom2, AtomPtr atom3)
{
	AtomType type1 = atom1->getGeomType();
	AtomType type2 = atom2->getGeomType();
	AtomType type3 = atom3->getGeomType();

	GeomTable table = GeomTable::getGeomTable();
	double angle = table.getBondAngle(type1, type2, type3);

	return angle;
}

void Atom::addProperties()
{
	addStringProperty("atom_name", &_atomName);
	addDoubleProperty("init_b", &_initialB);
	addVec3Property("init_pos", &_initialPosition);
	addVec3Property("pdb_pos", &_pdbPosition);
	addVec3Property("longest_axis", &_ellipsoidLongestAxis);
	addMat3x3Property("tensor", &_tensor);
	addIntProperty("atom_num", &_atomNum);
	addDoubleProperty("init_occupancy", &_origOccupancy);
	addStringProperty("conformer", &_conformer);
	addBoolProperty("from_pdb", &_fromPDB);
	addIntProperty("hetatm", &_hetatm);
	addDoubleProperty("weighting", &_weighting);

	if (_element)
	{
		_elementSymbol = _element->getSymbol();
	}

	addStringProperty("element", &_elementSymbol); 

	// add tensor, matrix stuff

	addChild("model", _model);
	addChild("dist_model", _distModelOnly);
}

void Atom::addObject(ParserPtr object, std::string category)
{
	if (category == "model")
	{
		ModelPtr model = ToModelPtr(object);
		setModel(model);
	}
	else if (category == "dist_model")
	{
		ModelPtr model = ToModelPtr(object);
		_distModelOnly = model;
	}
}

void Atom::postParseTidy()
{
	_element = Element::getElement(_elementSymbol);

	if (_element == ElementPtr())
	{
		std::cout << "Warning: element not found for " << shortDesc() << std::endl;
		std::cout << "Element symbol is: " << _elementSymbol << std::endl;
	}
	
	if (_hetatm < 0)
	{
		if (getModel()->isAbsolute())
		{
			_hetatm = ToAbsolutePtr(getModel())->isHeteroAtom();
		}
	}
}

bool Atom::closeToAtom(AtomPtr another, double tolerance)
{
	vec3 pos1 = getAbsolutePosition();
	vec3 pos2 = another->getAbsolutePosition();

	bool closeish = vec3_near_vec3_box(pos1, pos2, tolerance);

	if (!closeish)
	{
		return false;
	}
	else
	{
		vec3 diff = vec3_subtract_vec3(pos1, pos2);
		double length = vec3_length(diff);
		return (length < tolerance);
	}
}
