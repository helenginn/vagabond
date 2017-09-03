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

	addIdentityToType("asp", "C", AtomC);
	addIdentityToType("asp", "O", AtomO);
	addIdentityToType("asp", "N", AtomNH1);
	addIdentityToType("asp", "CA", AtomAspCA);
	addIdentityToType("asp", "CB", AtomAspCB);
	addIdentityToType("asp", "CG", AtomAspCG);
	addIdentityToType("asp", "OD1", AtomAspOD1);
	addIdentityToType("asp", "OD2", AtomAspOD2);

	addIdentityToType("leu", "C", AtomC);
	addIdentityToType("leu", "O", AtomO);
	addIdentityToType("leu", "N", AtomNH1);
	addIdentityToType("leu", "CA", AtomLeuCA);
	addIdentityToType("leu", "CB", AtomLeuCB);
	addIdentityToType("leu", "CG", AtomLeuCG);
	addIdentityToType("leu", "CD1", AtomLeuCD1);
	addIdentityToType("leu", "CD2", AtomLeuCD2);

	addIdentityToType("thr", "C", AtomC);
	addIdentityToType("thr", "O", AtomO);
	addIdentityToType("thr", "N", AtomNH1);
	addIdentityToType("thr", "CA", AtomThrCA);
	addIdentityToType("thr", "CB", AtomThrCB);
	addIdentityToType("thr", "OG1", AtomThrOG1);
	addIdentityToType("thr", "CG2", AtomThrCG2);

	addIdentityToType("cys", "C", AtomC);
	addIdentityToType("cys", "O", AtomO);
	addIdentityToType("cys", "N", AtomNH1);
	addIdentityToType("cys", "CA", AtomCysCA);
	addIdentityToType("cys", "CB", AtomCysCB);
	addIdentityToType("cys", "SG", AtomCysSG);

	addIdentityToType("ser", "C", AtomC);
	addIdentityToType("ser", "O", AtomO);
	addIdentityToType("ser", "N", AtomNH1);
	addIdentityToType("ser", "CA", AtomSerCA);
	addIdentityToType("ser", "CB", AtomSerCB);
	addIdentityToType("ser", "OG", AtomSerOG);

	addIdentityToType("val", "C", AtomC);
	addIdentityToType("val", "O", AtomO);
	addIdentityToType("val", "N", AtomNH1);
	addIdentityToType("val", "CA", AtomValCA);
	addIdentityToType("val", "CB", AtomValCB);
	addIdentityToType("val", "CG1", AtomValCG1);
	addIdentityToType("val", "CG2", AtomValCG2);

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
	addIdentityToType("his", "CA", AtomHisCA);
	addIdentityToType("his", "CB", AtomHisCB);
	addIdentityToType("his", "CG", AtomHisCG);
	addIdentityToType("his", "CE1", AtomHisCE1);
	addIdentityToType("his", "CD2", AtomHisCD2);
	addIdentityToType("his", "ND1", AtomHisND1);
	addIdentityToType("his", "NE2", AtomHisNE2);

	addIdentityToType("met", "C", AtomC);
	addIdentityToType("met", "O", AtomO);
	addIdentityToType("met", "N", AtomNH1);
	addIdentityToType("met", "CA", AtomMetCA);
	addIdentityToType("met", "CB", AtomMetCB);
	addIdentityToType("met", "CG", AtomMetCG);
	addIdentityToType("met", "SD", AtomMetSD);
	addIdentityToType("met", "CE", AtomMetCE);

	addIdentityToType("ile", "C", AtomC);
	addIdentityToType("ile", "O", AtomO);
	addIdentityToType("ile", "N", AtomNH1);
	addIdentityToType("ile", "CA", AtomIleCA);
	addIdentityToType("ile", "CB", AtomIleCB);
	addIdentityToType("ile", "CG1", AtomIleCG1);
	addIdentityToType("ile", "CG2", AtomIleCG2);
	addIdentityToType("ile", "CD1", AtomIleCD1);

	addIdentityToType("ala", "C", AtomC);
	addIdentityToType("ala", "O", AtomO);
	addIdentityToType("ala", "N", AtomNH1);
	addIdentityToType("ala", "CA", AtomAlaCA);
	addIdentityToType("ala", "CB", AtomAlaCB);

	addIdentityToType("gly", "C", AtomC);
	addIdentityToType("gly", "O", AtomO);
	addIdentityToType("gly", "N", AtomNH1);
	addIdentityToType("gly", "CA", AtomGlyCA);

	addIdentityToType("phe", "C", AtomC);
	addIdentityToType("phe", "O", AtomO);
	addIdentityToType("phe", "N", AtomNH1);
	addIdentityToType("phe", "CA", AtomPheCA);
	addIdentityToType("phe", "CB", AtomPheCB);
	addIdentityToType("phe", "CG", AtomPheCG);
	addIdentityToType("phe", "CD1", AtomPheCD1);
	addIdentityToType("phe", "CD2", AtomPheCD2);
	addIdentityToType("phe", "CE1", AtomPheCE1);
	addIdentityToType("phe", "CE2", AtomPheCE2);
	addIdentityToType("phe", "CZ", AtomPheCZ);

	addIdentityToType("tyr", "C", AtomC);
	addIdentityToType("tyr", "O", AtomO);
	addIdentityToType("tyr", "N", AtomNH1);
	addIdentityToType("tyr", "CA", AtomTyrCA);
	addIdentityToType("tyr", "CB", AtomTyrCB);
	addIdentityToType("tyr", "CG", AtomTyrCG);
	addIdentityToType("tyr", "CD1", AtomTyrCD1);
	addIdentityToType("tyr", "CD2", AtomTyrCD2);
	addIdentityToType("tyr", "CE1", AtomTyrCE1);
	addIdentityToType("tyr", "CE2", AtomTyrCE2);
	addIdentityToType("tyr", "CZ", AtomTyrCZ);
	addIdentityToType("tyr", "OH", AtomTyrOH);

	addIdentityToType("trp", "C", AtomC);
	addIdentityToType("trp", "O", AtomO);
	addIdentityToType("trp", "N", AtomNH1);
	addIdentityToType("trp", "CA", AtomTrpCA);
	addIdentityToType("trp", "CB", AtomTrpCB);
	addIdentityToType("trp", "CG", AtomTrpCG);
	addIdentityToType("trp", "CD1", AtomTrpCD1);
	addIdentityToType("trp", "CD2", AtomTrpCD2);
	addIdentityToType("trp", "NE1", AtomTrpNE1);
	addIdentityToType("trp", "CE2", AtomTrpCE2);
	addIdentityToType("trp", "CE3", AtomTrpCE3);
	addIdentityToType("trp", "CZ2", AtomTrpCZ2);
	addIdentityToType("trp", "CZ3", AtomTrpCZ3);
	addIdentityToType("trp", "CH2", AtomTrpCH2);

	addIdentityToType("pro", "C", AtomC);
	addIdentityToType("pro", "O", AtomO);
	addIdentityToType("pro", "N", AtomNH1);
	addIdentityToType("pro", "CA", AtomProCA);
	addIdentityToType("pro", "CB", AtomProCB);
	addIdentityToType("pro", "CG", AtomProCG);
	addIdentityToType("pro", "CD", AtomProCD);

	addIdentityToType("arg", "C", AtomC);
	addIdentityToType("arg", "O", AtomO);
	addIdentityToType("arg", "N", AtomNH1);
	addIdentityToType("arg", "CA", AtomArgCA);
	addIdentityToType("arg", "CB", AtomArgCB);
	addIdentityToType("arg", "CG", AtomArgCG);
	addIdentityToType("arg", "CD", AtomArgCD);
	addIdentityToType("arg", "NE", AtomArgNE);
	addIdentityToType("arg", "CZ", AtomArgCZ);
	addIdentityToType("arg", "NH1", AtomArgNH1);
	addIdentityToType("arg", "NH2", AtomArgNH2);

	addIdentityToType("asn", "C", AtomC);
	addIdentityToType("asn", "O", AtomO);
	addIdentityToType("asn", "N", AtomNH1);
	addIdentityToType("asn", "CA", AtomAsnCA);
	addIdentityToType("asn", "CB", AtomAsnCB);
	addIdentityToType("asn", "CG", AtomAsnCG);
	addIdentityToType("asn", "OD1", AtomAsnOD1);
	addIdentityToType("asn", "ND2", AtomAsnND2);

	addIdentityToType("gln", "C", AtomC);
	addIdentityToType("gln", "O", AtomO);
	addIdentityToType("gln", "N", AtomNH1);
	addIdentityToType("gln", "CA", AtomGlnCA);
	addIdentityToType("gln", "CB", AtomGlnCB);
	addIdentityToType("gln", "CG", AtomGlnCG);
	addIdentityToType("gln", "CD", AtomGlnCD);
	addIdentityToType("gln", "OE1", AtomGlnOE1);
	addIdentityToType("gln", "NE2", AtomGlnNE2);

	addIdentityToType("glu", "C", AtomC);
	addIdentityToType("glu", "O", AtomO);
	addIdentityToType("glu", "N", AtomNH1);
	addIdentityToType("glu", "CA", AtomGluCA);
	addIdentityToType("glu", "CB", AtomGluCB);
	addIdentityToType("glu", "CG", AtomGluCG);
	addIdentityToType("glu", "CD", AtomGluCD);
	addIdentityToType("glu", "OE1", AtomGluOE1);
	addIdentityToType("glu", "OE2", AtomGluOE2);

	addBondLength(AtomC, AtomGluCA, 1.525);
	addBondLength(AtomNH1, AtomGluCA, 1.459);
	addBondLength(AtomGluCA, AtomGluCB, 1.530);
	addBondLength(AtomGluCB, AtomGluCG, 1.520);
	addBondLength(AtomGluCG, AtomGluCD, 1.516);
	addBondLength(AtomGluCD, AtomGluOE1, 1.249);
	addBondLength(AtomGluCD, AtomGluOE2, 1.249);

	addBondAngle(AtomNH1, AtomGluCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomGluCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomGluCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomGluCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomGluCA, AtomGluCB, 110.5);
	addBondAngle(AtomGluCB, AtomGluCA, AtomC, 110.1);
	addBondAngle(AtomGluCA, AtomGluCB, AtomGluCG, 114.1);
	addBondAngle(AtomGluCB, AtomGluCG, AtomGluCD, 112.6);
	addBondAngle(AtomGluCG, AtomGluCD, AtomGluOE1, 118.4);
	addBondAngle(AtomGluCG, AtomGluCD, AtomGluOE2, 118.4);

	addBondLength(AtomC, AtomGlnCA, 1.525);
	addBondLength(AtomNH1, AtomGlnCA, 1.459);
	addBondLength(AtomGlnCA, AtomGlnCB, 1.530);
	addBondLength(AtomGlnCB, AtomGlnCG, 1.520);
	addBondLength(AtomGlnCG, AtomGlnCD, 1.516);
	addBondLength(AtomGlnCD, AtomGlnOE1, 1.231);
	addBondLength(AtomGlnCD, AtomGlnNE2, 1.328);

	addBondAngle(AtomNH1, AtomGlnCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomGlnCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomGlnCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomGlnCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomGlnCA, AtomGlnCB, 110.5);
	addBondAngle(AtomGlnCB, AtomGlnCA, AtomC, 110.1);
	addBondAngle(AtomGlnCA, AtomGlnCB, AtomGlnCG, 114.1);
	addBondAngle(AtomGlnCB, AtomGlnCG, AtomGlnCD, 112.6);
	addBondAngle(AtomGlnCG, AtomGlnCD, AtomGlnOE1, 120.8);
	addBondAngle(AtomGlnCG, AtomGlnCD, AtomGlnNE2, 116.4);

	addBondLength(AtomC, AtomAsnCA, 1.525);
	addBondLength(AtomNH1, AtomAsnCA, 1.459);
	addBondLength(AtomAsnCA, AtomAsnCB, 1.530);
	addBondLength(AtomAsnCB, AtomAsnCG, 1.516);
	addBondLength(AtomAsnCG, AtomAsnOD1, 1.231);
	addBondLength(AtomAsnCG, AtomAsnND2, 1.328);

	addBondAngle(AtomNH1, AtomAsnCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomAsnCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomAsnCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomAsnCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomAsnCA, AtomAsnCB, 110.5);
	addBondAngle(AtomAsnCB, AtomAsnCA, AtomC, 110.1);
	addBondAngle(AtomAsnCA, AtomAsnCB, AtomAsnCG, 112.6);
	addBondAngle(AtomAsnCB, AtomAsnCG, AtomAsnOD1, 120.8);
	addBondAngle(AtomAsnCB, AtomAsnCG, AtomAsnND2, 116.4);

	addBondLength(AtomC, AtomArgCA, 1.525);
	addBondLength(AtomNH1, AtomArgCA, 1.459);
	addBondLength(AtomArgCA, AtomArgCB, 1.530);
	addBondLength(AtomArgCB, AtomArgCG, 1.520);
	addBondLength(AtomArgCG, AtomArgCD, 1.520);
	addBondLength(AtomArgCD, AtomArgNE, 1.460);
	addBondLength(AtomArgNE, AtomArgCZ, 1.329);
	addBondLength(AtomArgCZ, AtomArgNH1, 1.326);
	addBondLength(AtomArgCZ, AtomArgNH2, 1.326);

	addBondAngle(AtomNH1, AtomArgCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomArgCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomArgCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomArgCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomArgCA, AtomArgCB, 110.5);
	addBondAngle(AtomArgCB, AtomArgCA, AtomC, 110.1);
	addBondAngle(AtomArgCA, AtomArgCB, AtomArgCG, 114.1);
	addBondAngle(AtomArgCB, AtomArgCG, AtomArgCD, 111.3);
	addBondAngle(AtomArgCG, AtomArgCD, AtomArgNE, 112.0);
	addBondAngle(AtomArgCD, AtomArgNE, AtomArgCZ, 124.2);
	addBondAngle(AtomArgNE, AtomArgCZ, AtomArgNH1, 120.0);
	addBondAngle(AtomArgNE, AtomArgCZ, AtomArgNH2, 119.7);

	addBondAngle(AtomNH1, AtomProCA, AtomC, 112.1);
	addBondAngle(AtomProCA, AtomC, AtomO, 120.2);
	addBondAngle(AtomC, AtomNH1, AtomProCA, 119.3);
	addBondAngle(AtomProCA, AtomC, AtomNH1, 117.1);
	addBondAngle(AtomNH1, AtomProCA, AtomProCB, 103.0);
	addBondAngle(AtomProCB, AtomProCA, AtomC, 110.1);
	addBondAngle(AtomProCA, AtomProCB, AtomProCG, 104.5);
	addBondAngle(AtomProCB, AtomProCG, AtomProCD, 106.1);

	addBondLength(AtomC, AtomProCA, 1.524);
	addBondLength(AtomNH1, AtomProCA, 1.468);
	addBondLength(AtomProCA, AtomProCB, 1.530);
	addBondLength(AtomProCB, AtomProCG, 1.492);
	addBondLength(AtomProCG, AtomProCD, 1.503);
	
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
	addBondAngle(AtomMetCB, AtomMetCG, AtomMetSD, 112.7);
	addBondAngle(AtomMetCG, AtomMetSD, AtomMetCE, 100.9);

	addBondLength(AtomC, AtomGlyCA, 1.514);
	addBondLength(AtomNH1, AtomGlyCA, 1.456);

	addBondAngle(AtomNH1, AtomGlyCA, AtomC, 113.1);
	addBondAngle(AtomGlyCA, AtomC, AtomNH1, 116.2);
	addBondAngle(AtomGlyCA, AtomC, AtomO, 120.6);
	addBondAngle(AtomGlyCA, AtomNH1, AtomC, 122.3);

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
	addBondAngle(AtomNH1, AtomLysCA, AtomLysCB, 110.5);
	addBondAngle(AtomLysCB, AtomLysCA, AtomC, 110.1);
	addBondAngle(AtomLysCA, AtomLysCB, AtomLysCG, 114.1);
	addBondAngle(AtomLysCB, AtomLysCG, AtomLysCD, 111.3);
	addBondAngle(AtomLysCG, AtomLysCD, AtomLysCE, 111.3);
	addBondAngle(AtomLysCD, AtomLysCE, AtomLysNZ, 111.9);

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

	addBondLength(AtomC, AtomCysCA, 1.525);
	addBondLength(AtomNH1, AtomCysCA, 1.459);
	addBondLength(AtomCysCA, AtomCysCB, 1.530);
	addBondLength(AtomCysCB, AtomCysSG, 1.808);

	addBondAngle(AtomNH1, AtomCysCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomCysCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomCysCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomCysCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomCysCA, AtomCysCB, 110.5);
	addBondAngle(AtomCysCB, AtomCysCA, AtomC, 110.1);
	addBondAngle(AtomCysCA, AtomCysCB, AtomCysSG, 114.4);
	
	addBondLength(AtomC, AtomHisCA, 1.525);
	addBondLength(AtomNH1, AtomHisCA, 1.459);
	addBondLength(AtomHisCA, AtomHisCB, 1.530);
	addBondLength(AtomHisCB, AtomHisCG, 1.497);
	addBondLength(AtomHisCG, AtomHisND1, 1.371);
	addBondLength(AtomHisCG, AtomHisCD2, 1.356);
	addBondLength(AtomHisND1, AtomHisCE1, 1.319);
	addBondLength(AtomHisCD2, AtomHisNE2, 1.374);

	addBondAngle(AtomNH1, AtomHisCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomHisCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomHisCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomHisCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomHisCA, AtomHisCB, 110.5);
	addBondAngle(AtomHisCB, AtomHisCA, AtomC, 110.1);
	addBondAngle(AtomHisCA, AtomHisCB, AtomHisCG, 113.8);
	addBondAngle(AtomHisCB, AtomHisCG, AtomHisCD2, 131.2);
	addBondAngle(AtomHisCB, AtomHisCG, AtomHisND1, 122.7);
	addBondAngle(AtomHisCG, AtomHisCD2, AtomHisNE2, 107.2);
	addBondAngle(AtomHisCG, AtomHisND1, AtomHisCE1, 109.3);

	addBondLength(AtomC, AtomPheCA, 1.525);
	addBondLength(AtomNH1, AtomPheCA, 1.459);
	addBondLength(AtomPheCA, AtomPheCB, 1.530);
	addBondLength(AtomPheCB, AtomPheCG, 1.502);
	addBondLength(AtomPheCG, AtomPheCD1, 1.384);
	addBondLength(AtomPheCG, AtomPheCD2, 1.384);
	addBondLength(AtomPheCD1, AtomPheCE1, 1.384);
	addBondLength(AtomPheCD2, AtomPheCE2, 1.384);
	addBondLength(AtomPheCE2, AtomPheCZ, 1.382);

	addBondAngle(AtomNH1, AtomPheCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomPheCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomPheCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomPheCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomPheCA, AtomPheCB, 110.5);
	addBondAngle(AtomPheCB, AtomPheCA, AtomC, 110.1);
	addBondAngle(AtomPheCA, AtomPheCB, AtomPheCG, 113.8);
	addBondAngle(AtomPheCB, AtomPheCG, AtomPheCD2, 120.7);
	addBondAngle(AtomPheCB, AtomPheCG, AtomPheCD1, 120.7);
	addBondAngle(AtomPheCG, AtomPheCD2, AtomPheCE2, 120.7);
	addBondAngle(AtomPheCG, AtomPheCD1, AtomPheCE1, 120.7);
	addBondAngle(AtomPheCD2, AtomPheCE2, AtomPheCZ, 120.0);

	addBondAngle(AtomNH1, AtomIleCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomIleCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomIleCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomIleCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomIleCA, AtomIleCB, 111.5);
	addBondAngle(AtomIleCB, AtomIleCA, AtomC, 109.1);
	addBondAngle(AtomIleCA, AtomIleCB, AtomIleCG1, 110.4);
	addBondAngle(AtomIleCA, AtomIleCB, AtomIleCG2, 110.5);
	addBondAngle(AtomIleCB, AtomIleCG1, AtomIleCD1, 113.8);

	addBondLength(AtomC, AtomIleCA, 1.525);
	addBondLength(AtomNH1, AtomIleCA, 1.459);
	addBondLength(AtomIleCA, AtomIleCB, 1.540);
	addBondLength(AtomIleCB, AtomIleCG1, 1.530);
	addBondLength(AtomIleCG1, AtomIleCD1, 1.513);
	addBondLength(AtomIleCB, AtomIleCG2, 1.521);

	addBondLength(AtomC, AtomAlaCA, 1.525);
	addBondLength(AtomNH1, AtomAlaCA, 1.459);
	addBondLength(AtomAlaCA, AtomAlaCB, 1.521);

	addBondAngle(AtomNH1, AtomAlaCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomAlaCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomAlaCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomAlaCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomAlaCA, AtomAlaCB, 110.4);
	addBondAngle(AtomAlaCB, AtomAlaCA, AtomC, 110.5);

	addBondLength(AtomC, AtomLeuCA, 1.525);
	addBondLength(AtomNH1, AtomLeuCA, 1.459);
	addBondLength(AtomLeuCA, AtomLeuCB, 1.530);
	addBondLength(AtomLeuCB, AtomLeuCG, 1.530);
	addBondLength(AtomLeuCG, AtomLeuCD1, 1.521);
	addBondLength(AtomLeuCG, AtomLeuCD2, 1.521);

	addBondAngle(AtomNH1, AtomLeuCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomLeuCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomLeuCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomLeuCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomLeuCA, AtomLeuCB, 110.5);
	addBondAngle(AtomLeuCB, AtomLeuCA, AtomC, 110.1);
	addBondAngle(AtomLeuCA, AtomLeuCB, AtomLeuCG, 116.3);
	addBondAngle(AtomLeuCB, AtomLeuCG, AtomLeuCD1, 110.7);
	addBondAngle(AtomLeuCB, AtomLeuCG, AtomLeuCD2, 110.7);

	addBondLength(AtomC, AtomTrpCA, 1.525);
	addBondLength(AtomNH1, AtomTrpCA, 1.459);
	addBondLength(AtomTrpCA, AtomTrpCB, 1.530);
	addBondLength(AtomTrpCB, AtomTrpCG, 1.498);
	addBondLength(AtomTrpCG, AtomTrpCD1, 1.365);
	addBondLength(AtomTrpCG, AtomTrpCD2, 1.433);
	addBondLength(AtomTrpCD1, AtomTrpNE1, 1.374);
	addBondLength(AtomTrpCD2, AtomTrpCE2, 1.409);
	addBondLength(AtomTrpCD2, AtomTrpCE3, 1.398);
	addBondLength(AtomTrpCE2, AtomTrpCZ2, 1.394);
	addBondLength(AtomTrpCE3, AtomTrpCZ3, 1.382);
	addBondLength(AtomTrpCZ3, AtomTrpCH2, 1.368);
	addBondLength(AtomTrpCZ2, AtomTrpCH2, 1.400);

	addBondAngle(AtomNH1, AtomTrpCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomTrpCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomTrpCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomTrpCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomTrpCA, AtomTrpCB, 110.5);
	addBondAngle(AtomTrpCB, AtomTrpCA, AtomC, 110.1);
	addBondAngle(AtomTrpCA, AtomTrpCB, AtomTrpCG, 113.6);
	addBondAngle(AtomTrpCB, AtomTrpCG, AtomTrpCD1, 126.9);
	addBondAngle(AtomTrpCB, AtomTrpCG, AtomTrpCD2, 126.8);
	addBondAngle(AtomTrpCG, AtomTrpCD1, AtomTrpNE1, 110.2);
	addBondAngle(AtomTrpCG, AtomTrpCD2, AtomTrpCE3, 133.9);
	addBondAngle(AtomTrpCD1, AtomTrpNE1, AtomTrpCE2, 108.9);
	addBondAngle(AtomTrpCZ3, AtomTrpCE3, AtomTrpCD2, 118.6);
	addBondAngle(AtomTrpCD2, AtomTrpCE2, AtomTrpCZ2, 122.4);

	addBondLength(AtomC, AtomTyrCA, 1.525);
	addBondLength(AtomNH1, AtomTyrCA, 1.459);
	addBondLength(AtomTyrCA, AtomTyrCB, 1.530);
	addBondLength(AtomTyrCB, AtomTyrCG, 1.512);
	addBondLength(AtomTyrCG, AtomTyrCD1, 1.389);
	addBondLength(AtomTyrCG, AtomTyrCD2, 1.389);
	addBondLength(AtomTyrCD1, AtomTyrCE1, 1.382);
	addBondLength(AtomTyrCD2, AtomTyrCE2, 1.382);
	addBondLength(AtomTyrCE2, AtomTyrCZ, 1.378);
	addBondLength(AtomTyrCZ, AtomTyrOH, 1.376);

	addBondAngle(AtomNH1, AtomTyrCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomTyrCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomTyrCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomTyrCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomTyrCA, AtomTyrCB, 110.5);
	addBondAngle(AtomTyrCB, AtomTyrCA, AtomC, 110.1);
	addBondAngle(AtomTyrCA, AtomTyrCB, AtomTyrCG, 113.9);
	addBondAngle(AtomTyrCB, AtomTyrCG, AtomTyrCD2, 120.8);
	addBondAngle(AtomTyrCB, AtomTyrCG, AtomTyrCD1, 120.8);
	addBondAngle(AtomTyrCG, AtomTyrCD2, AtomTyrCE2, 121.2);
	addBondAngle(AtomTyrCG, AtomTyrCD1, AtomTyrCE1, 121.2);
	addBondAngle(AtomTyrCD2, AtomTyrCE2, AtomTyrCZ, 119.6);
	addBondAngle(AtomTyrCE2, AtomTyrCZ, AtomTyrOH, 119.9);

	addBondLength(AtomC, AtomAspCA, 1.525);
	addBondLength(AtomNH1, AtomAspCA, 1.459);
	addBondLength(AtomAspCA, AtomAspCB, 1.530);
	addBondLength(AtomAspCB, AtomAspCG, 1.516);
	addBondLength(AtomAspCG, AtomAspOD1, 1.249);
	addBondLength(AtomAspCG, AtomAspOD2, 1.249);

	addBondAngle(AtomNH1, AtomAspCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomAspCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomAspCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomAspCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomAspCA, AtomAspCB, 110.5);
	addBondAngle(AtomAspCB, AtomAspCA, AtomC, 110.1);
	addBondAngle(AtomAspCA, AtomAspCB, AtomAspCG, 112.6);
	addBondAngle(AtomAspCB, AtomAspCG, AtomAspOD1, 118.4);
	addBondAngle(AtomAspCB, AtomAspCG, AtomAspOD2, 118.4);

	addBondAngle(AtomNH1, AtomAspCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomAspCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomAspCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomAspCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomAspCA, AtomAspCB, 110.5);
	addBondAngle(AtomAspCB, AtomAspCA, AtomC, 110.1);
	addBondAngle(AtomAspCA, AtomAspCB, AtomAspCG, 112.6);
	addBondAngle(AtomAspCB, AtomAspCG, AtomAspOD1, 118.4);
	addBondAngle(AtomAspCB, AtomAspCG, AtomAspOD2, 118.4);

	addBondLength(AtomC, AtomThrCA, 1.525);
	addBondLength(AtomNH1, AtomThrCA, 1.459);
	addBondLength(AtomThrCA, AtomThrCB, 1.540);
	addBondLength(AtomThrCB, AtomThrOG1, 1.433);
	addBondLength(AtomThrCB, AtomThrCG2, 1.521);

	addBondAngle(AtomNH1, AtomThrCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomThrCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomThrCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomThrCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomThrCA, AtomThrCB, 111.5);
	addBondAngle(AtomThrCB, AtomThrCA, AtomC, 109.1);
	addBondAngle(AtomThrCA, AtomThrCB, AtomThrCG2, 110.5);
	addBondAngle(AtomThrCA, AtomThrCB, AtomThrOG1, 109.6);

	addBondLength(AtomC, AtomValCA, 1.525);
	addBondLength(AtomNH1, AtomValCA, 1.459);
	addBondLength(AtomValCA, AtomValCB, 1.540);
	addBondLength(AtomValCB, AtomValCG1, 1.524);
	addBondLength(AtomValCB, AtomValCG2, 1.524);

	addBondAngle(AtomNH1, AtomValCA, AtomC, C_CH1E_NH1_ANGLE);
	addBondAngle(AtomValCA, AtomC, AtomO, CH1E_C_O_ANGLE);
	addBondAngle(AtomC, AtomNH1, AtomValCA, C_NH1_CH1E_ANGLE);
	addBondAngle(AtomValCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
	addBondAngle(AtomNH1, AtomValCA, AtomValCB, 111.5);
	addBondAngle(AtomValCB, AtomValCA, AtomC, 109.1);
	addBondAngle(AtomValCA, AtomValCB, AtomValCG1, 110.5);
	addBondAngle(AtomValCA, AtomValCB, AtomValCG2, 110.5);

}

