//
//  GeomTable.cpp
//  vagabond
//
//  Created by Helen Ginn on 07/08/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

// (1) Acta Crystallographica Section D: Biological Crystallography, 63(5), 611-620.

#define CH2E_CH2E_LENGTH (1.520)
#define CH2E_NH3_LENGTH  (1.489)
#define CH1E_OH1_LENGTH  (1.433)
#define CH2E_OH1_LENGTH  (1.417)
#define CH1E_CH2E_LENGTH (1.530)
#define CH1E_CH3E_LENGTH (1.521)
#define CH1E_CH1E_LENGTH (1.540)
#define C_CH1E_LENGTH    (1.523) // // 1.525 - Engh&Huber, 1.523 - ref1.
#define CH1E_NH1_LENGTH  (1.455) // 1.458 - Engh&Huber, 1.455 - ref1.
#define C_NH1_LENGTH     (1.332) // 1.329 - E&H - used ref1.
#define C_O_LENGTH       (1.231) // agrees with ref1.
#define CH2E_SM_LENGTH   (1.803)
#define CH3E_SM_LENGTH   (1.791)
#define C5_CH2E_LENGTH   (1.497)
#define C5_CR1H_LENGTH   (1.354)
#define C5_NH1_LENGTH    (1.378)
#define CRH_NH1_LENGTH   (1.345)
#define CR1H_NH1_LENGTH  (1.374)
#define CR1E_NH1_LENGTH  (1.374)

#define CH1E_CH2E_CH2E_ANGLE 114.1 
#define CH2E_CH2E_CH2E_ANGLE 111.3 
#define CH1E_CH1E_CH3E_ANGLE 110.5

#define CH2E_CH2E_NH3_ANGLE  111.9
#define CH1E_CH1E_OH1_ANGLE  109.6 
#define CH1E_CH2E_OH1_ANGLE  111.1
#define CH2E_C5_CR1H_ANGLE   131.2
#define CH2E_C5_CR1E_ANGLE   129.1
#define CH2E_C5_NH1_ANGLE    122.7
#define CR1H_C5_NH1_ANGLE    106.1
#define CR1E_C5_NH1_ANGLE    105.2
#define C5_CH2E_CH1E_ANGLE   113.8
#define C5_CR1H_NH1_ANGLE    107.2
#define CH2E_SM_CH3E_ANGLE   100.2
#define CH2E_CH2E_SM_ANGLE   112.4// looks wrong
#define CH1E_C_NH1_ANGLE     117.2
#define CH1E_C_O_ANGLE       120.1
#define NH1_C_O_ANGLE        122.7
#define C_CH1E_NH1_ANGLE     111.0 // ref1 says too high! 111.2 - E&H
#define NH1_CH1E_CH2E_ANGLE  110.6
#define C_NH1_CH1E_ANGLE     121.7

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
	addBondLength(AtomCH2E, AtomOH1, CH2E_OH1_LENGTH);
	addBondLength(AtomCH1E, AtomCH2E, CH1E_CH2E_LENGTH);
	addBondLength(AtomCH1E, AtomCH3E, CH1E_CH3E_LENGTH);
	addBondLength(AtomCH1E, AtomCH1E, CH1E_CH1E_LENGTH);
	addBondLength(AtomC5, AtomCH2E, C5_CH2E_LENGTH);
	addBondLength(AtomC5, AtomCR1H, C5_CR1H_LENGTH);
	addBondLength(AtomC5, AtomNH1, C5_NH1_LENGTH);
	addBondLength(AtomCRH, AtomNH1, CRH_NH1_LENGTH);
	addBondLength(AtomCR1H, AtomNH1, CR1H_NH1_LENGTH);
	addBondLength(AtomCR1E, AtomNH1, CR1E_NH1_LENGTH);
	addBondLength(AtomCH1E, AtomC, C_CH1E_LENGTH);
	addBondLength(AtomCH1E, AtomNH1, CH1E_NH1_LENGTH);
	addBondLength(AtomNH1, AtomC, C_NH1_LENGTH);
	addBondLength(AtomC, AtomO, C_O_LENGTH);
	addBondLength(AtomCH2E, AtomSM, CH2E_SM_LENGTH);
	addBondLength(AtomCH3E, AtomSM, CH3E_SM_LENGTH);

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
	addBondAngle(AtomCH2E, AtomSM, AtomCH3E, CH2E_SM_CH3E_ANGLE);
	addBondAngle(AtomCH2E, AtomCH2E, AtomSM, CH2E_CH2E_SM_ANGLE);

	addBondAngle(AtomCH1E, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomCH1E, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomNH1, AtomC, AtomO, NH1_C_O_ANGLE);
	addBondAngle(AtomC, AtomCH1E, AtomNH1, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomCH2E, AtomCH1E, AtomNH1, NH1_CH1E_CH2E_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomCH1E, C_NH1_CH1E_ANGLE);


	addIdentityToType("thr", "C", AtomC);
	addIdentityToType("thr", "O", AtomO);
	addIdentityToType("thr", "N", AtomNH1);
	addIdentityToType("thr", "CA", AtomCH1E);
	addIdentityToType("thr", "CB", AtomCH1E);
	addIdentityToType("thr", "OG1", AtomOH1);
	addIdentityToType("thr", "CG2", AtomCH3E);

	addIdentityToType("ser", "C", AtomC);
	addIdentityToType("ser", "O", AtomO);
	addIdentityToType("ser", "N", AtomNH1);
	addIdentityToType("ser", "CA", AtomSerCA);
	addIdentityToType("ser", "CB", AtomSerCB);
	addIdentityToType("ser", "OG", AtomSerOG);

	addIdentityToType("val", "C", AtomC);
	addIdentityToType("val", "O", AtomO);
	addIdentityToType("val", "N", AtomNH1);
	addIdentityToType("val", "CA", AtomCH1E);
	addIdentityToType("val", "CB", AtomCH1E);
	addIdentityToType("val", "CG1", AtomCH3E);
	addIdentityToType("val", "CG2", AtomCH3E);

	addIdentityToType("lys", "C", AtomC);
	addIdentityToType("lys", "O", AtomO);
	addIdentityToType("lys", "N", AtomNH1);
	addIdentityToType("lys", "CA", AtomLysCA);
	addIdentityToType("lys", "CB", AtomLysCB);
	addIdentityToType("lys", "CG", AtomLysCG);
	addIdentityToType("lys", "CD", AtomLysCD);
	addIdentityToType("lys", "CE", AtomLysCE);
	addIdentityToType("lys", "NZ", AtomLysNZ);

	addIdentityToType("his", "C", AtomC);
	addIdentityToType("his", "O", AtomO);
	addIdentityToType("his", "N", AtomNH1);
	addIdentityToType("his", "CA", AtomCH1E);
	addIdentityToType("his", "CB", AtomCH2E);
	addIdentityToType("his", "CG", AtomC5);
	addIdentityToType("his", "CE1", AtomCR1E);
	addIdentityToType("his", "CD2", AtomCR1H);
	addIdentityToType("his", "ND1", AtomNH1);
	addIdentityToType("his", "NE2", AtomNH1);

	addIdentityToType("met", "C", AtomC);
	addIdentityToType("met", "O", AtomO);
	addIdentityToType("met", "N", AtomNH1);
	addIdentityToType("met", "CA", AtomMetCA);
	addIdentityToType("met", "CB", AtomMetCB);
	addIdentityToType("met", "CG", AtomMetCG);
	addIdentityToType("met", "SD", AtomMetSD);
	addIdentityToType("met", "CE", AtomMetCE);

	addIdentityToType("phe", "C", AtomC);
	addIdentityToType("phe", "O", AtomO);
	addIdentityToType("phe", "N", AtomNH1);
	addIdentityToType("phe", "CA", AtomCH1E);
	addIdentityToType("phe", "CB", AtomCH2E);

	addIdentityToType("ile", "N", AtomNH1);

	addBondLength(AtomC, AtomMetCA, 1.525);
	addBondLength(AtomNH1, AtomMetCA, 1.459);
	addBondLength(AtomMetCA, AtomMetCB, 1.535);
	addBondLength(AtomMetCB, AtomMetCG, 1.509);
	addBondLength(AtomMetCG, AtomMetSD, 1.807);
	addBondLength(AtomMetSD, AtomMetCE, 1.774);


	addBondAngle(AtomNH1, AtomMetCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomMetCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomMetCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomMetCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomMetCA, AtomMetCB, 110.5);
	addBondAngle(AtomMetCA, AtomMetCB, AtomC, 110.1);
	addBondAngle(AtomMetCA, AtomMetCB, AtomMetCG, 114.1);
	addBondAngle(AtomMetCB, AtomMetCG, AtomMetSD, 116.0);
	addBondAngle(AtomMetCG, AtomMetSD, AtomMetCE, 100.9);

	addBondLength(AtomC, AtomLysCA, 1.525);
	addBondLength(AtomNH1, AtomLysCA, 1.459);
	addBondLength(AtomLysCA, AtomLysCB, 1.535);
	addBondLength(AtomLysCB, AtomLysCG, 1.521);
	addBondLength(AtomLysCG, AtomLysCD, 1.520);
	addBondLength(AtomLysCD, AtomLysCE, 1.508);
	addBondLength(AtomLysCE, AtomLysNZ, 1.486);

	addBondAngle(AtomNH1, AtomLysCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomLysCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomLysCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomLysCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomLysCA, AtomLysCB, 110.6);
	addBondAngle(AtomLysCB, AtomLysCA, AtomC, 110.4);
	addBondAngle(AtomLysCA, AtomLysCB, AtomLysCG, 113.4);
	addBondAngle(AtomLysCB, AtomLysCG, AtomLysCD, 111.6);
	addBondAngle(AtomLysCG, AtomLysCD, AtomLysCE, 111.9);
	addBondAngle(AtomLysCD, AtomLysCE, AtomLysNZ, 111.7);

	addBondLength(AtomC, AtomSerCA, 1.525);
	addBondLength(AtomNH1, AtomSerCA, 1.459);
	addBondLength(AtomSerCA, AtomSerCB, 1.525);
	addBondLength(AtomSerCB, AtomSerOG, 1.418);

	addBondAngle(AtomNH1, AtomSerCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomSerCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomSerCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomSerCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomSerCA, AtomSerCB, 110.5);
	addBondAngle(AtomSerCB, AtomSerCA, AtomC, 110.1);
	addBondAngle(AtomSerCA, AtomSerCB, AtomSerOG, 111.2);

}

