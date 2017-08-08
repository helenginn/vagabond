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

#define CH1E_CH2E_CH2E_ANGLE 114.1 // 2.0 sigma
#define CH2E_CH2E_CH2E_ANGLE 111.3 // 2.3 sigma
#define CH2E_CH2E_NH3_ANGLE  111.9 // 3.2 sigma
#define CH1E_CH1E_OH1_ANGLE  109.6 // 1.5 sigma
#define CH1E_CH1E_CH3E_ANGLE 110.5 // 3.0 sigma

#include "shared_ptrs.h"
#include "GeomTable.h"

GeomTable GeomTable::_geomTable;

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


GeomTable::GeomTable()
{
	addBondLength(AtomCH2E, AtomCH2E, CH2E_CH2E_LENGTH);
	addBondLength(AtomCH2E, AtomNH3, CH2E_NH3_LENGTH);
	addBondLength(AtomCH1E, AtomOH1, CH1E_OH1_LENGTH);
	addBondLength(AtomCH1E, AtomCH2E, CH1E_CH2E_LENGTH);
	addBondLength(AtomCH1E, AtomCH3E, CH1E_CH3E_LENGTH);
	addBondLength(AtomCH1E, AtomCH1E, CH1E_CH1E_LENGTH);

	addBondAngle(AtomCH1E, AtomCH2E, AtomCH2E, CH1E_CH2E_CH2E_ANGLE);
	addBondAngle(AtomCH2E, AtomCH2E, AtomCH2E, CH2E_CH2E_CH2E_ANGLE);
	addBondAngle(AtomCH2E, AtomCH2E, AtomNH3, CH2E_CH2E_NH3_ANGLE);
	addBondAngle(AtomCH1E, AtomCH1E, AtomOH1, CH1E_CH1E_OH1_ANGLE);
	addBondAngle(AtomCH1E, AtomCH1E, AtomCH3E, CH1E_CH1E_CH3E_ANGLE);

}