//
//  GeomTable.cpp
//  vagabond
//
//  Created by Helen Ginn on 07/08/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#define CH2E_CH2E_LENGTH (1.520)
#define CH2E_NH3_LENGTH  (1.489)
#define CH1E_OH1_LENGTH  (1.433)
#define CH1E_CH2E_LENGTH (1.530)
#define CH1E_CH3E_LENGTH (1.521)
#define CH1E_CH1E_LENGTH (1.540)

#define C5_CH2E_LENGTH (1.497)
#define C5_CR1H_LENGTH (1.354)
#define C5_NH1_LENGTH (1.378)
#define CRH_NH1_LENGTH (1.345)
#define CR1H_NH1_LENGTH (1.374)
#define CR1E_NH1_LENGTH (1.374)

#define CH1E_CH2E_CH2E_ANGLE 114.1 
#define CH2E_CH2E_CH2E_ANGLE 111.3 
#define CH2E_CH2E_NH3_ANGLE  111.9 
#define CH1E_CH1E_OH1_ANGLE  109.6 
#define CH1E_CH2E_OH1_ANGLE  111.1
#define CH1E_CH1E_CH3E_ANGLE 110.5
#define CH2E_C5_CR1H_ANGLE   131.2 
#define CH2E_C5_CR1E_ANGLE   129.1
#define CH2E_C5_NH1_ANGLE    122.7
#define CR1H_C5_NH1_ANGLE    106.1
#define CR1E_C5_NH1_ANGLE    105.2
#define C5_CH2E_CH1E_ANGLE   113.8
#define C5_CR1H_NH1_ANGLE    107.2

#include "shared_ptrs.h"
#include "GeomTable.h"

GeomTable GeomTable::_geomTable;

AtomType GeomTable::getType(std::string res, std::string atomName)
{
	AtomIdentity ident;
	ident.first = res;
	ident.second = atomName;

	if (!_identityToType.count(ident))
	{
		return AtomUnassigned;
	}

	return _identityToType[ident];
}

double GeomTable::getBondLength(AtomType atom1, AtomType atom2)
{
	AtomPair pair;
	pair.first = atom1;
	pair.second = atom2;

	if (_bondLengths.count(pair))
	{
		return _bondLengths[pair];
	}

	return -1;
}

double GeomTable::getBondAngle(AtomType atom1, AtomType atom2, AtomType atom3)
{
	AtomPair pair;
	pair.first = atom1;
	pair.second = atom2;

	AtomTrio trio;
	trio.first = pair;
	trio.second = atom3;

	if (_bondAngles.count(trio))
	{
		return _bondAngles[trio];
	}

	return -1;
}

void GeomTable::addBondLength(AtomType atom1, AtomType atom2, double length)
{
	AtomPair pair;
	pair.first = atom1;
	pair.second = atom2;

	_bondLengths[pair] = length;

	AtomPair reverse;
	reverse.second = atom1;
	reverse.first = atom2;

	_bondLengths[reverse] = length;
}

void GeomTable::addBondAngle(AtomType atom1, AtomType atom2,
							 AtomType atom3, double angle)
{
	AtomTrio trio;
	trio.first.first = atom1;
	trio.first.second = atom2;
	trio.second = atom3;

	_bondAngles[trio] = deg2rad(angle);

	AtomTrio reverse;
	reverse.first.first = atom3;
	reverse.first.second = atom2;
	reverse.second = atom1;

	_bondAngles[reverse] = deg2rad(angle);
}

void GeomTable::addIdentityToType(std::string res, std::string atomName,
								  AtomType type)
{
	AtomIdentity ident;
	ident.first = res;
	ident.second = atomName;

	_identityToType[ident] = type;
}

GeomTable::GeomTable()
{
	addBondLength(AtomCH2E, AtomCH2E, CH2E_CH2E_LENGTH);
	addBondLength(AtomCH2E, AtomNH3, CH2E_NH3_LENGTH);
	addBondLength(AtomCH1E, AtomOH1, CH1E_OH1_LENGTH);
	addBondLength(AtomCH1E, AtomCH2E, CH1E_CH2E_LENGTH);
	addBondLength(AtomCH1E, AtomCH3E, CH1E_CH3E_LENGTH);
	addBondLength(AtomCH1E, AtomCH1E, CH1E_CH1E_LENGTH);
	addBondLength(AtomC5, AtomCH2E, C5_CH2E_LENGTH);
	addBondLength(AtomC5, AtomCR1H, C5_CR1H_LENGTH);
	addBondLength(AtomC5, AtomNH1, C5_NH1_LENGTH);
	addBondLength(AtomCRH, AtomNH1, CRH_NH1_LENGTH);
	addBondLength(AtomCR1H, AtomNH1, CR1H_NH1_LENGTH);
	addBondLength(AtomCR1E, AtomNH1, CR1E_NH1_LENGTH);

	addBondAngle(AtomCH1E, AtomCH2E, AtomCH2E, CH1E_CH2E_CH2E_ANGLE);
	addBondAngle(AtomCH2E, AtomCH2E, AtomCH2E, CH2E_CH2E_CH2E_ANGLE);
	addBondAngle(AtomCH2E, AtomCH2E, AtomNH3, CH2E_CH2E_NH3_ANGLE);
	addBondAngle(AtomCH1E, AtomCH1E, AtomOH1, CH1E_CH1E_OH1_ANGLE);
	addBondAngle(AtomCH1E, AtomCH1E, AtomCH3E, CH1E_CH1E_CH3E_ANGLE);
	addBondAngle(AtomCH1E, AtomCH2E, AtomOH1, CH1E_CH2E_OH1_ANGLE);
	addBondAngle(AtomCH2E, AtomC5, AtomCR1H, CH2E_C5_CR1H_ANGLE);
	addBondAngle(AtomCH2E, AtomC5, AtomCR1E, CH2E_C5_CR1E_ANGLE);
	addBondAngle(AtomCH2E, AtomC5, AtomNH1, CH2E_C5_NH1_ANGLE);
	addBondAngle(AtomCR1H, AtomC5, AtomNH1, CR1H_C5_NH1_ANGLE);
	addBondAngle(AtomCR1E, AtomC5, AtomNH1, CR1E_C5_NH1_ANGLE);
	addBondAngle(AtomC5, AtomCH2E, AtomCH1E, C5_CH2E_CH1E_ANGLE);
	addBondAngle(AtomC5, AtomCR1H, AtomNH1, C5_CR1H_NH1_ANGLE);

	addIdentityToType("thr", "CA", AtomCH1E);
	addIdentityToType("thr", "CB", AtomCH1E);
	addIdentityToType("thr", "OG1", AtomOH1);
	addIdentityToType("thr", "CG2", AtomCH3E);

	addIdentityToType("ser", "CA", AtomCH1E);
	addIdentityToType("ser", "CB", AtomCH2E);
	addIdentityToType("ser", "OG", AtomOH1);

	addIdentityToType("val", "CA", AtomCH1E);
	addIdentityToType("val", "CB", AtomCH1E);
	addIdentityToType("val", "CG1", AtomCH3E);
	addIdentityToType("val", "CG2", AtomCH3E);

	addIdentityToType("lys", "CA", AtomCH1E);
	addIdentityToType("lys", "CB", AtomCH2E);
	addIdentityToType("lys", "CG", AtomCH2E);
	addIdentityToType("lys", "CD", AtomCH2E);
	addIdentityToType("lys", "CE", AtomCH2E);
	addIdentityToType("lys", "NZ", AtomNH3);

	addIdentityToType("his", "CA", AtomCH1E);
	addIdentityToType("his", "CB", AtomCH2E);
	addIdentityToType("his", "CG", AtomC5);
	addIdentityToType("his", "CE1", AtomCR1E);
	addIdentityToType("his", "CD2", AtomCR1H);
	addIdentityToType("his", "ND1", AtomNH1);
	addIdentityToType("his", "NE2", AtomNH1);
}

