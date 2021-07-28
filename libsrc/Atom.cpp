//
//  Atom.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Atom.h"
#include "Absolute.h"
#include "FFT.h"
#include <hcsrc/mat3x3.h>
#include <cmath>
#include <stdlib.h>
#include "Bond.h"
#include "Whack.h"
#include <iostream>
#include "Element.h"
#include "Monomer.h"
#include "Polymer.h"
#include <iomanip>
#include <hcsrc/FileReader.h>
#include <sstream>
#include "Shouter.h"
#include "Crystal.h"
#include "PDBReader.h"
#include "Anisotropicator.h"
#include "GhostBond.h"
#include "WaterNetwork.h"

Atom::Atom()
{
	_asu = -1;
	_initialPosition = make_vec3(0, 0, 0);
	_pdbPosition = make_vec3(0, 0, 0);
	_ellipsoidLongestAxis = make_vec3(0, 0, 0);
	_initialB = 0;
	_geomType = AtomUnassigned;
	_weighting = 1;
	_targetWeight = 0;
	_targetPos = empty_vec3();
	_origOccupancy = 1.0;
	_fromPDB = true;
	_isWater = 0;
	_tensor = make_mat3x3();
	_hetatm = 0;
	_hBondage = false;
	_weightOnly = 1;
}

Atom::Atom(Atom &other)
{
	_tensor = other._tensor;
	_initialPosition = other._initialPosition;
	_initialB = other._initialB;
	_geomType = other._geomType;
	_element = other._element;
	_atomName = other._atomName;
	_model = other._model;
	_monomer = other._monomer;
	_weighting = 1;
	_weightOnly = 1;
	_atomNum = other._atomNum;
	_fromPDB = other._fromPDB;
	_pdbPosition = other._pdbPosition;
	_origOccupancy = other._origOccupancy;
	_conformer = other._conformer;
	_ellipsoidLongestAxis = other._ellipsoidLongestAxis;
	_hetatm = other._hetatm;
}

Atom::~Atom()
{

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
		if (!getMonomer())
		{
			return AtomUnassigned;
		}
		std::string id = getMonomer()->getIdentifier();
		findAtomType(id);
	}

	return _geomType;
}

double Atom::getSolventRadius()
{
	AtomType geom = getGeomType();

	// carbonyl carbons and ring carbons are 2.1
	if (geom == AtomC || geom == AtomProC ||
		geom == AtomTyrCD1 || geom == AtomTyrCD2 ||
		geom == AtomTyrCE1 || geom == AtomTyrCE2 ||
		geom == AtomTyrCZ ||
		geom == AtomPheCD1 || geom == AtomPheCD2 ||
		geom == AtomPheCE1 || geom == AtomPheCE2 ||
		geom == AtomHisCE1 || geom == AtomHisCD2)
	{
		return 2.1;
	}
	else if (getElement()->getSymbol() == "C")
	{
		return 2.3;
	}
	else if (getElement()->getSymbol() == "N")
	{
		return 1.6;
	}
	else if (getElement()->getSymbol() == "O")
	{
		return 1.6;
	}
	else if (getElement()->getSymbol() == "S")
	{
		return 1.9;
	}
	else if (getElement()->getSymbol() == "H")
	{
		return 0;
	}
	
	return 2.0;
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
	
	refreshBondAngles();
}

void Atom::refreshBondAngles()
{
	BondPtr bond = ToBondPtr(getModel());
	ModelPtr parent = bond->getParentModel();
	
	if (!parent->isBond())
	{
		return;	
	}
	
	ToBondPtr(parent)->resetBondAngles();

}

std::string Atom::description()
{
	std::ostringstream ss;
	ss << std::endl;
	ss << "Atom " << shortDesc() << std::endl;
	ss << "  Model type: " << getModel()->getClassName() << std::endl;
	ss << "  From PDB? " << (_fromPDB ? "yes" : "no") << std::endl;
	
	if (_fromPDB)
	{
		ss << "  PDB initial position: " << 
		vec3_desc(_initialPosition) << std::endl;
	}
	
	ss << "  Is HETATM? " << (_hetatm ? "yes" : "no") << std::endl;
	ss << std::endl;
	
	return ss.str();
}

void Atom::inheritParents()
{
	getMonomer()->addAtom(shared_from_this());
	getMonomer()->getPolymer()->addAtom(shared_from_this());
}

bool Atom::isWater()
{
	return (getMolecule() && getMolecule()->isWaterNetwork());
}

void Atom::setModel(ModelPtr model)
{
	if (_model)
	{
		model->setMolecule(_model->getMolecule());
	}

	_model = model;
}

vec3 Atom::getSymRelatedPosition(int i, vec3 pos)
{
	CrystalPtr crystal = DynCrystalPtr(getTopParser());
	CSym::CCP4SPG *spg = crystal->getSpaceGroup();
	
	mat3x3 real2Frac = crystal->getReal2Frac();
	mat3x3_mult_vec(real2Frac, &pos);

	float *rot = &spg->symop[i].rot[0][0];
	float *trn = spg->symop[i].trn;

	vec3 mod = empty_vec3();
	mod.x = pos.x * rot[0] + pos.y * rot[1] + pos.z * rot[2];
	mod.y = pos.x * rot[3] + pos.y * rot[4] + pos.z * rot[5];
	mod.z = pos.x * rot[6] + pos.y * rot[7] + pos.z * rot[8];
	mod.x += trn[0];
	mod.y += trn[1];
	mod.z += trn[2];
	
	fmod(mod.x, 1);
	fmod(mod.y, 1);
	fmod(mod.z, 1);
	
	mat3x3 frac2Real = crystal->getFrac2Real();
	mat3x3_mult_vec(frac2Real, &mod);
	
	return mod;
}

vec3 Atom::getSymRelatedPosition(int i, int conf)
{
	vec3 pos = getAbsolutePosition();
	
	if (conf >= 0 && getModel()->hasExplicitPositions())
	{
		pos = getExplicitModel()->getFinalPositions()[conf].start;
	}
	
	return getSymRelatedPosition(i, pos);
}

vec3 Atom::getPositionInUnitCell()
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	CSym::CCP4SPG *spg = crystal->getSpaceGroup();
	mat3x3 f2r = crystal->getFFT()->toRecip();
	mat3x3 r2f = crystal->getFFT()->toReal();

	vec3 tmp = getAbsolutePosition();
	mat3x3_mult_vec(f2r, &tmp);
	VagFFT::collapseFrac(&tmp.x, &tmp.y, &tmp.z);
	mat3x3_mult_vec(r2f, &tmp);

	return tmp;
}

vec3 Atom::getPositionInAsu(int conf)
{
	CrystalPtr crystal = DynCrystalPtr(getTopParser());
	CSym::CCP4SPG *spg = crystal->getSpaceGroup();
	mat3x3 f2rt = crystal->getReal2Frac();
	mat3x3 f2r = mat3x3_transpose(f2rt);
	mat3x3 r2f = mat3x3_inverse(f2r);
	int nsymop = spg->nsymop;
	
	for (int i = 0; i < nsymop; i++)
	{
		vec3 pos = getSymRelatedPosition(i, conf);
		vec3 tmp = pos;
		mat3x3_mult_vec(f2r, &tmp);
		VagFFT::collapseFrac(&tmp.x, &tmp.y, &tmp.z);
		
		if (tmp.x < spg->mapasu_zero[0] &&
		    tmp.y < spg->mapasu_zero[1] &&
		    tmp.z < spg->mapasu_zero[2]) 
		{
			mat3x3_mult_vec(r2f, &tmp);
			return tmp;
		}
	}
	
	return empty_vec3();
}

size_t Atom::symOpCount()
{
	CrystalPtr crystal = DynCrystalPtr(getTopParser());
	CSym::CCP4SPG *spg = crystal->getSpaceGroup();
	return spg->nsymop;
}

void Atom::addPointerToLocalArea(VagFFTPtr fft, vec3 pos,
                                 std::vector<Atom *> *ptrs, double rad)
{
	/* Limiting case must be far enough away to merge 'nicely' into the
	 * constant-1 solvent mask area */

//	double b = this->getBFactor();
	mat3x3 unit_cell = fft->getRecipBasis();
	mat3x3_mult_vec(unit_cell, &pos);
	
	vec3 mins, maxs;
	mat3x3 basis = fft->getRealBasis();	
	fft->findLimitingValues(-rad, rad, -rad, 
	                        rad, -rad, rad,
	                        &mins, &maxs);
	
	for (int k = mins.z; k < maxs.z; k++)
	{
		for (int j = mins.y; j < maxs.y; j++)
		{
			for (int i = mins.x; i < maxs.x; i++)
			{
				long ele = fft->element(pos.x+i, pos.y+j, pos.z+k);
				Atom *current = ptrs->at(ele);
				
//				if (current && current->getBFactor() > b)
				{
//					continue;
				}
				
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(basis, &ijk);
				
				if (vec3_sqlength(ijk) > rad * rad)
				{
					continue;
				}
				
				ptrs->at(ele) = this;
			}
		}
	}
}

void Atom::addManyToMask(VagFFTPtr fft, int conf, int total)
{

	if (getElectronCount() <= 1 || _weighting <= 0)
	{
		return;
	}

	double radius = getSolventRadius();
	
	if (radius <= 0)
	{
		return;
	}
	
	radius += Options::getActiveCrystal()->getProbeRadius();

	for (int i = 0; i < total; i++)
	{
		if (getModel()->hasExplicitPositions() && 
		    i + conf > getExplicitModel()->getFinalPositions().size())
		{
			break;
		}

		size_t max = symOpCount();
		for (int j = 0; j < max; j++)
		{
			vec3 pos = empty_vec3();

			if (!getModel()->hasExplicitPositions())
			{
				pos = getSymRelatedPosition(j, -1);
			}
			else
			{
				pos = getSymRelatedPosition(j, i + conf);
			}

			fft->addToValueAroundPoint(pos, radius, 1, i);
		}
	}
}

void Atom::addToSolventMask(VagFFTPtr fft, double rad,
							std::vector<Atom *> *ptrs, int conf)
{
	if (getElectronCount() <= 1 || _weighting <= 0)
	{
		return;
	}

	double radius = getSolventRadius();
	
	if (radius <= 0)
	{
		return;
	}
	
	radius += Options::getActiveCrystal()->getProbeRadius();

	vec3 pos = getAbsolutePosition();
	
	if (conf >= 0 && getModel()->hasExplicitPositions())
	{
		pos = getExplicitModel()->getFinalPositions()[conf].start;
	}
	
	if (conf < 0)
	{
		size_t max = symOpCount();
		for (size_t i = 0; i < max; i++)
		{
			vec3 pos = getSymRelatedPosition(i, -1);
			fft->addToValueAroundPoint(pos, radius, -1);
	
			addPointerToLocalArea(fft, pos, ptrs, radius);
		}
	}
}

vec3 Atom::getAbsolutePosition()
{
	return _model->getAbsolutePosition();
}

ExplicitModelPtr Atom::getExplicitModel()
{
	return ToExplicitModelPtr(getModel());
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
	_geomType = GeomTable::getGeomTable()->getType(resName, _atomName);
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

	if (isHeteroAtom() && start == "ATOM  ")
	{
		start = "HETATM";
		resNum = _atomNum;
	}


	if (getMolecule())
	{
		chainID = getMolecule()->getChainID()[0];
	}
	
	if (getModel()->isAbsolute())
	{
		chainID = ToAbsolutePtr(getModel())->getChainID()[0];
		resNum = ToAbsolutePtr(getModel())->getResNum();
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
	
	if (_atomName.length() <= 3)
	{
		line << "  " << std::left << std::setfill(' ') << std::setw(3) << _atomName;
	}
	else
	{
		line << " " << std::left << std::setfill(' ') << std::setw(4) << _atomName;
	}

	line << std::right << conformer;
	line << std::setw(3) << residueName;
	line << " " << chainID;
	line << std::setfill(' ') << std::setw(4) << resNum;
	line << "  ";

	return line.str();
}

double Atom::posToMouse()
{
	getModel()->refreshPositions();
	vec3 bestPos = getModel()->getAbsolutePosition();
	vec3 target = _targetPos;

	vec3 diff = vec3_subtract_vec3(bestPos, target);
	double score = vec3_length(diff) * _targetWeight;

	return score;

}

void Atom::saveInitialPosition()
{
	vec3 pos = getModel()->getAbsolutePosition();
	_initialPosition = pos;
}

double Atom::posDisplacement(bool fromSaved, bool refresh, bool sq)
{
	if (!isFromPDB() && !fromSaved)
	{
		return 0;
	}
	
	if (fromSaved && getElectronCount() == 1)
	{
		return 0;
	}

	if (refresh)
	{
		getModel()->refreshPositions();
	}

	vec3 bestPos = getModel()->getAbsolutePosition();
	vec3 initialPos = getPDBPosition();
	
	if (fromSaved)
	{
		initialPos = getInitialPosition();
	}

	vec3 diff = vec3_subtract_vec3(bestPos, initialPos);
	double score = vec3_length(diff);
	
	if (sq)
	{
		score *= score;
	}

	return score;
}

std::string Atom::anisouPDBLine(CrystalPtr)
{
	ElementPtr element = getElement();
	if (element->getSymbol() == "H")
	{
		return "";
	}

	std::ostringstream stream;
	stream << pdbLineBeginning("ANISOU");

	mat3x3 realTensor = getModel()->getRealSpaceTensor();

	double scale = 10e4;
	stream << std::setprecision(0) << std::fixed;
	stream << std::setfill(' ') << std::setw(7) << scale * realTensor.vals[0];
	stream << std::setfill(' ') << std::setw(7) << scale * realTensor.vals[4];
	stream << std::setfill(' ') << std::setw(7) << scale * realTensor.vals[8];
	stream << std::setfill(' ') << std::setw(7) << scale * realTensor.vals[1];
	stream << std::setfill(' ') << std::setw(7) << scale * realTensor.vals[2];
	stream << std::setfill(' ') << std::setw(7) << scale * realTensor.vals[5];
	stream << std::endl;

	return stream.str();
}

std::string Atom::averagePDBContribution(bool samePos, bool sameB)
{
	getModel()->refreshPositions();
	std::string atomName = getAtomName();

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

	if (!getModel()->hasExplicitPositions())
	{
		return "";
	}
	
	std::ostringstream stream;
	std::vector<BondSample> samples = getExplicitModel()->getFinalPositions();
	std::vector<vec3> positions = getExplicitModel()->fishPositions(NULL);

	int i = ensembleNum;
	int count = 0;

	vec3 placement = positions[i];
	double occupancy = samples[i].occupancy * positions.size();
	stream << PDBReader::writeLine(shared_from_this(), placement, 
	                               count, occupancy, bFac);

	return stream.str();
}

std::string Atom::longDesc()
{
	std::string str = getMolecule()->getChainID() + "_" + shortDesc();

	return str;
}

std::string Atom::shortDesc()
{
	if (!getMonomer())
	{
		return getAtomName() + "_" + i_to_str(getAtomNum()) + _conformer;
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
	
	if (getModel() && getModel()->getMolecule())
	{
		return getModel()->getMolecule();
	}

	return MoleculePtr();
}

double Atom::getAngle(AtomPtr atom1, AtomPtr atom2, AtomPtr atom3)
{
	AtomType type1 = atom1->getGeomType();
	AtomType type2 = atom2->getGeomType();
	AtomType type3 = atom3->getGeomType();

	GeomTable *table = GeomTable::getGeomTable();
	double angle = table->getBondAngle(type1, type2, type3);

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
	addBoolProperty("water", &_isWater);
	addDoubleProperty("weighting", &_weighting);

	if (_element)
	{
		_elementSymbol = _element->getSymbol();
	}

	addStringProperty("element", &_elementSymbol); 

	// add tensor, matrix stuff

	addChild("model", _model);
	addChild("ghost", _ghost);
	addChild("dist_model", _distModelOnly);
}

void Atom::addObject(ParserPtr object, std::string category)
{
	if (category == "model")
	{
		ModelPtr model = ToModelPtr(object);
		setModel(model);
	}
	else if (category == "ghost")
	{
		GhostBondPtr ghost = ToGhostBondPtr(object);
		setGhostBond(ghost);
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

	if (!getModel()->hasExplicitPositions())
	{
		return -1;
	}

	ExplicitModelPtr expModel = ToExplicitModelPtr(getModel());
	ExplicitModelPtr otherExp = ToExplicitModelPtr(other->getModel());

	if (nSample >= 0)
	{
		me = expModel->getFinalPositions()[nSample].start;
		you = otherExp->getFinalPositions()[nSample].start;
	}
	
	if (!quick)
	{
		me = getAbsolutePosition();
		you = other->getAbsolutePosition();
	}
	
	vec3 apart = vec3_subtract_vec3(me, you);
	
	double dist = vec3_length(apart);
	
	return dist;
}

void Atom::writePositionsToFile(std::string suffix)
{
	std::string filename = "atom_positions_" + shortDesc() + suffix + ".csv";
	
	getExplicitModel()->writePositionsToFile(filename);
}

int Atom::getElectronCount()
{
	return _element->electronCount();
}

bool Atom::isAtom(std::string atomName, int resNum)
{
	return (getMonomer() && getMonomer()->getResidueNum() == resNum 
	        && _atomName == atomName);
}

double Atom::fishWhackMagnitude()
{
	std::lock_guard<std::mutex> lock(_whackLock);

	MonomerPtr m = getMonomer();

	if (!m)
	{
		return -1;
	}

	AtomPtr ca = m->findAtom("CA");

	if (!ca || !ca->getModel()->isBond())
	{
		return -1;
	}

	BondPtr b = ToBondPtr(ca->getModel());

	if (!b->hasWhack())
	{
		return -1;
	}

	WhackPtr whack = b->getWhack();
	double tot = fabs(Whack::getWhack(&*whack));
	tot += fabs(Whack::getKick(&*whack));

	return tot;

	// is sidechain

}

bool Atom::greater(AtomPtr a1, AtomPtr a2)
{
	double z1 = a1->getAbsolutePosition().z;
	double z2 = a2->getAbsolutePosition().z;

	if (a1->getExplicitModel()->hasExplicitPositions())
	{
		z1 = a1->getExplicitModel()->getLowestZ();
	}

	if (a2->getExplicitModel()->hasExplicitPositions())
	{
		z2 = a2->getExplicitModel()->getLowestZ();
	}
	
	return (z2 > z1);
}

std::string Atom::getChainID()
{
	if (!getMolecule())
	{
		return "";
	}
	
	return getMolecule()->getChainID();
}

void Atom::setElement(ElementPtr element)
{
	_element = element;
	_elementSymbol = element->getSymbol();
}
