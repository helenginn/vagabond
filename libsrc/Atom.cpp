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
#include "Plucker.h"
#include <iomanip>
#include "FileReader.h"
#include <sstream>
#include "Crystal.h"
#include "PDBReader.h"
#include "Anisotropicator.h"

Atom::Atom()
{
	_asu = -1;
	_waterPlucker = NULL;
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
	_hBondage = false;
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

void Atom::addToMap(FFTPtr fft, mat3x3 unit_cell, vec3 offset, bool mask,
                    bool sameScale, bool noWrap)
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

	MapScoreType type = (noWrap ? MapScoreAddNoWrap : MapScoreTypeNone);
	
	int solvent = Options::getAddSolvent();

	if (!mask || (mask && solvent < 2))
	{
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

		FFT::operation(fft, modified, pos, type, NULL, sameScale);
	}
	else
	{
		std::vector<BondSample> samples = _model->getFinalPositions();

		/* Each addition should only contribute the right occupancy */
		modified->multiplyAll(1 / (double)samples.size());
		
		/* Loop through each position and add separately */
		
		for (int i = 0; i < samples.size(); i++)
		{
			vec3 pos = samples[i].start;
			pos = vec3_subtract_vec3(pos, offset);
			mat3x3_mult_vec(unit_cell, &pos);

			FFT::operation(fft, modified, pos, type, NULL, sameScale);
		}
	}
}

vec3 Atom::getAbsolutePosition()
{
	return _model->getAbsolutePosition();
}

/* Need to add symops. */
vec3 Atom::getAsymUnitPosition(CrystalPtr crystal, int nSample)
{
	vec3 pos = getAbsolutePosition();
	
	if (nSample >= 0 && nSample < getModel()->getFinalPositions().size())
	{
		pos = getModel()->getFinalPositions()[nSample].start;
	}
	
	mat3x3 real2frac = crystal->getReal2Frac();
	mat3x3 frac2real = crystal->getHKL2Real();

	mat3x3_mult_vec(real2frac, &pos);
	
	if (pos.x != pos.x || !isfinite(pos.x))
	{
		return pos;
	}
	
	FFT::collapseFrac(&pos.x, &pos.y, &pos.z);
	CSym::CCP4SPG *spg = crystal->getSpaceGroup();
	pos = FFT::collapseToRealASU(pos, spg);
	mat3x3_mult_vec(frac2real, &pos);
	
	return pos;
}

bool Atom::isBackbone()
{
	if (_atomName == "C") return true;
	if (_atomName == "H") return true;
	if (_atomName == "N") return true;
	if (_atomName == "O") return true;
	if (_atomName == "OXT") return true;
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
	std::string residueName = "HOH";

	int resNum = 1;
	std::string chainID = "Z";
	
	if (getMonomer())
	{
		residueName = getMonomer()->getIdentifier();
		resNum = getMonomer()->getResidueNum();
	}

	if (isHeteroAtom())
	{
		start = "HETATM";
		resNum = _atomNum;
	}


	if (getMolecule())
	{
		chainID = getMolecule()->getChainID()[0];
	}

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
	line << " " << chainID;
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

double Atom::getBFactor()
{
	return _model->getMeanSquareDeviation();
}

std::string Atom::getPDBContribution(int ensembleNum)
{
	std::string atomName = getAtomName();
	ElementPtr element = getElement();
	double bFac = Options::getRuntimeOptions()->getGlobalBFactor();

	int tries = 10;

	if (element->getSymbol() == "H")
	{
		return "";
		tries = 1;
	}

	std::ostringstream stream;
	std::vector<BondSample> positions = getModel()->getFinalPositions();

	int i = ensembleNum;
	int count = 0;

	vec3 placement = positions[i].start;
	double occupancy = positions[i].occupancy * positions.size();
	stream << PDBReader::writeLine(shared_from_this(), placement, 
	                               count, occupancy, bFac);

	return stream.str();
}

void Atom::cacheCloseWaters(double tolerance)
{
	if (_waterPlucker != NULL)
	{
		delete _waterPlucker;
		_waterPlucker = NULL;
	}
	
	_waterPlucker = new Plucker();
	_waterPlucker->setGranularity(0.2);
	
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	
	std::vector<AtomPtr> atoms = crystal->getCloseAtoms(shared_from_this(),
	                                                    tolerance);
	
	for (int i = 0; i < atoms.size(); i++)
	{
		AtomPtr atm = atoms[i];
		
		if (atm->getAtomName() != "O")
		{
			continue;
		}
		
		double occ = atm->getModel()->getEffectiveOccupancy();
		_waterPlucker->addPluckable(&*atm, occ);
	}
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
	addBoolProperty("hbonding", &_hBondage);
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

double Atom::getDistanceFrom(Atom *other, int nSample, bool quick)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	vec3 me = getAbsolutePosition();
	vec3 you = other->getAbsolutePosition();

	if (nSample >= 0)
	{
		me = getModel()->getFinalPositions()[nSample].start;
		you = other->getModel()->getFinalPositions()[nSample].start;
	}
	
	if (!quick)
	{
		me = getAsymUnitPosition(crystal, nSample);
		you = other->getAsymUnitPosition(crystal, nSample);
	}
	
	vec3 apart = vec3_subtract_vec3(me, you);
	
	double dist = vec3_length(apart);
	
	return dist;
}

Atom *Atom::pluckAnother()
{
	if (_waterPlucker == NULL)
	{
		return NULL;
	}
	
	return static_cast<Atom *>(_waterPlucker->pluck());
}

size_t Atom::pluckCount()
{
	if (_waterPlucker == NULL)
	{
		return 0;
	}

	return _waterPlucker->pluckCount();
}
