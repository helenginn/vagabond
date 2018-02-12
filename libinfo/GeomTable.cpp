//
//  GeomTable.cpp
//  vagabond
//
//  Created by Helen Ginn on 07/08/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

// Including corrections from
// (1) Acta Crystallographica Section D: Biological Crystallography, 63(5), 611-620.
//

#define CH2E_CH2E_LENGTH (1.520)
#define CH2E_NH3_LENGTH  (1.489)
#define CH1E_OH1_LENGTH  (1.433)
#define CH2E_OH1_LENGTH  (1.417)
#define CH1E_CH2E_LENGTH (1.530)
#define CH1E_CH3E_LENGTH (1.521)
#define CH1E_CH1E_LENGTH (1.540)

#define C_CH1E_LENGTH    (1.523) // // 1.525 - Engh&Huber, 1.523 - ref1.
#define CH1E_NH1_LENGTH  (1.455) // 1.458 - Engh&Huber, 1.455 - ref1.
#define C_NH1_LENGTH     (1.329) // 1.329 - E&H, 1.332 - ref1.
#define C_O_LENGTH       (1.231) // 1.229 - E&H, 1.231 - ref1.

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

#define CH1E_C_NH1_ANGLE     117.2 // 117.2
#define CH1E_C_O_ANGLE       120.1 // 120.1 - Engh&Huber, 120.7 - CCP4
#define C_NH1_CH1E_ANGLE     121.7
#define NH1_C_O_ANGLE        122.7
#define C_CH1E_NH1_ANGLE     111.2 // matches Engh & Huber and CCP4

#define deg2rad(a) ((a) * M_PI / 180)

#include "GeomTable.h"
#include <math.h>
#include <iostream>

GeomTable GeomTable::_geomTable;

AtomType GeomTable::getType(std::string res, std::string atomName)
{
    AtomIdentity ident;
    ident.first = res;
    ident.second = atomName;

    if (!_identityToType.count(ident))
    {
        return AtomUnknown;
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
    _three2OneCode["pro"] = "P";
    _three2OneCode["cys"] = "C";
    _three2OneCode["met"] = "M";
    _three2OneCode["tyr"] = "Y";
    _three2OneCode["trp"] = "W";
    _three2OneCode["phe"] = "F";
    _three2OneCode["his"] = "H";
    _three2OneCode["ile"] = "I";
    _three2OneCode["leu"] = "L";
    _three2OneCode["val"] = "V";
    _three2OneCode["thr"] = "T";
    _three2OneCode["ser"] = "S";
    _three2OneCode["asp"] = "D";
    _three2OneCode["asn"] = "N";
    _three2OneCode["glu"] = "E";
    _three2OneCode["gln"] = "Q";
    _three2OneCode["lys"] = "K";
    _three2OneCode["arg"] = "R";
    _three2OneCode["ala"] = "A";
    _three2OneCode["gly"] = "G";

    addBondLength(AtomCH2E, AtomCH2E, CH2E_CH2E_LENGTH);
    addBondLength(AtomCH2E, AtomNH3, CH2E_NH3_LENGTH);
    addBondLength(AtomCH1E, AtomC, C_CH1E_LENGTH);
    addBondLength(AtomCH1E, AtomNH1, CH1E_NH1_LENGTH);
    addBondLength(AtomNH1, AtomC, C_NH1_LENGTH);
    addBondLength(AtomC, AtomO, C_O_LENGTH);
    addBondAngle(AtomCH1E, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomCH1E, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomNH1, AtomC, AtomO, NH1_C_O_ANGLE);
    addBondAngle(AtomC, AtomCH1E, AtomNH1, C_CH1E_NH1_ANGLE);
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
    addBondLength(AtomGluCA, AtomGluCB, 1.535);
    addBondLength(AtomGluCB, AtomGluCG, 1.517);
    addBondLength(AtomGluCG, AtomGluCD, 1.515);
    addBondLength(AtomGluCD, AtomGluOE1, 1.252);
    addBondLength(AtomGluCD, AtomGluOE2, 1.252);

    addBondAngle(AtomNH1, AtomGluCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomGluCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomGluCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomGluCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomGluCA, AtomGluCB, 110.6);
    addBondAngle(AtomGluCB, AtomGluCA, AtomC, 110.4);
    addBondAngle(AtomGluCA, AtomGluCB, AtomGluCG, 113.4);
    addBondAngle(AtomGluCB, AtomGluCG, AtomGluCD, 114.2);
    addBondAngle(AtomGluCG, AtomGluCD, AtomGluOE1, 118.3);
    addBondAngle(AtomGluCG, AtomGluCD, AtomGluOE2, 118.3);

    addBondLength(AtomC, AtomGlnCA, 1.525);
    addBondLength(AtomNH1, AtomGlnCA, 1.459);
    addBondLength(AtomGlnCA, AtomGlnCB, 1.535);
    addBondLength(AtomGlnCB, AtomGlnCG, 1.521);
    addBondLength(AtomGlnCG, AtomGlnCD, 1.506);
    addBondLength(AtomGlnCD, AtomGlnOE1, 1.235);
    addBondLength(AtomGlnCD, AtomGlnNE2, 1.324);

    addBondAngle(AtomNH1, AtomGlnCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomGlnCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomGlnCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomGlnCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomGlnCA, AtomGlnCB, 110.6);
    addBondAngle(AtomGlnCB, AtomGlnCA, AtomC, 110.4);
    addBondAngle(AtomGlnCA, AtomGlnCB, AtomGlnCG, 113.4);
    addBondAngle(AtomGlnCB, AtomGlnCG, AtomGlnCD, 111.6);
    addBondAngle(AtomGlnCG, AtomGlnCD, AtomGlnOE1, 121.6);
    addBondAngle(AtomGlnCG, AtomGlnCD, AtomGlnNE2, 116.7);

    addBondLength(AtomC, AtomAsnCA, 1.525);
    addBondLength(AtomNH1, AtomAsnCA, 1.459);
    addBondLength(AtomAsnCA, AtomAsnCB, 1.527);
    addBondLength(AtomAsnCB, AtomAsnCG, 1.506);
    addBondLength(AtomAsnCG, AtomAsnOD1, 1.235);
    addBondLength(AtomAsnCG, AtomAsnND2, 1.324);

    addBondAngle(AtomNH1, AtomAsnCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomAsnCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomAsnCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomAsnCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomAsnCA, AtomAsnCB, 110.6);
    addBondAngle(AtomAsnCB, AtomAsnCA, AtomC, 110.4);
    addBondAngle(AtomAsnCA, AtomAsnCB, AtomAsnCG, 113.4);
    addBondAngle(AtomAsnCB, AtomAsnCG, AtomAsnOD1, 121.6);
    addBondAngle(AtomAsnCB, AtomAsnCG, AtomAsnND2, 116.6);
    addBondAngle(AtomAsnOD1, AtomAsnCG, AtomAsnND2, 121.8);

    addBondLength(AtomC, AtomArgCA, 1.525);
    addBondLength(AtomNH1, AtomArgCA, 1.459);
    addBondLength(AtomArgCA, AtomArgCB, 1.535);
    addBondLength(AtomArgCB, AtomArgCG, 1.521);
    addBondLength(AtomArgCG, AtomArgCD, 1.515);
    addBondLength(AtomArgCD, AtomArgNE, 1.460);
    addBondLength(AtomArgNE, AtomArgCZ, 1.326);
    addBondLength(AtomArgCZ, AtomArgNH1, 1.326);
    addBondLength(AtomArgCZ, AtomArgNH2, 1.326);

    addBondAngle(AtomNH1, AtomArgCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomArgCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomArgCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomArgCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomArgCA, AtomArgCB, 110.6);
    addBondAngle(AtomArgCB, AtomArgCA, AtomC, 110.4);
    addBondAngle(AtomArgCA, AtomArgCB, AtomArgCG, 113.4);
    addBondAngle(AtomArgCB, AtomArgCG, AtomArgCD, 111.6);
    addBondAngle(AtomArgCG, AtomArgCD, AtomArgNE, 111.8);
    addBondAngle(AtomArgCD, AtomArgNE, AtomArgCZ, 123.6);
    addBondAngle(AtomArgNE, AtomArgCZ, AtomArgNH1, 120.3);
    addBondAngle(AtomArgNE, AtomArgCZ, AtomArgNH2, 119.4);

    addBondAngle(AtomNH1, AtomProCA, AtomC, 112.48);
    addBondAngle(AtomProCA, AtomC, AtomO, 120.2);
    addBondAngle(AtomC, AtomNH1, AtomProCA, 119.3);
    addBondAngle(AtomProCA, AtomC, AtomNH1, 117.1);
    addBondAngle(AtomNH1, AtomProCA, AtomProCB, 103.0);
    addBondAngle(AtomProCB, AtomProCA, AtomC, 110.1);
    addBondAngle(AtomProCA, AtomProCB, AtomProCG, 104.5);
    addBondAngle(AtomProCB, AtomProCG, AtomProCD, 106.1);

    addBondLength(AtomC, AtomProCA, 1.524);
    addBondLength(AtomNH1, AtomProCA, 1.468);
    addBondLength(AtomProCA, AtomProCB, 1.531);
    addBondLength(AtomProCB, AtomProCG, 1.495);
    addBondLength(AtomProCG, AtomProCD, 1.502);

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
    addBondAngle(AtomNH1, AtomMetCA, AtomMetCB, 110.6);
    addBondAngle(AtomC, AtomMetCA, AtomMetCB, 110.4);
    addBondAngle(AtomMetCA, AtomMetCB, AtomMetCG, 113.3);
    addBondAngle(AtomMetCB, AtomMetCG, AtomMetSD, 112.4);
    addBondAngle(AtomMetCG, AtomMetSD, AtomMetCE, 100.2);

    addBondLength(AtomC, AtomGlyCA, 1.514);
    addBondLength(AtomNH1, AtomGlyCA, 1.456);

    addBondAngle(AtomNH1, AtomGlyCA, AtomC, 113.9);
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

    addBondLength(AtomC, AtomCysCA, 1.525);
    addBondLength(AtomNH1, AtomCysCA, 1.459);
    addBondLength(AtomCysCA, AtomCysCB, 1.526);
    addBondLength(AtomCysCB, AtomCysSG, 1.812);

    addBondAngle(AtomNH1, AtomCysCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomCysCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomCysCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomCysCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomCysCA, AtomCysCB, 110.8);
    addBondAngle(AtomCysCB, AtomCysCA, AtomC, 111.5);
    addBondAngle(AtomCysCA, AtomCysCB, AtomCysSG, 114.2);

    addBondLength(AtomC, AtomHisCA, 1.525);
    addBondLength(AtomNH1, AtomHisCA, 1.459);
    addBondLength(AtomHisCA, AtomHisCB, 1.535);
    addBondLength(AtomHisCB, AtomHisCG, 1.492);
    addBondLength(AtomHisCG, AtomHisND1, 1.380);
    addBondLength(AtomHisCG, AtomHisCD2, 1.354);
    addBondLength(AtomHisND1, AtomHisCE1, 1.326);
    addBondLength(AtomHisCD2, AtomHisNE2, 1.373);

    addBondAngle(AtomNH1, AtomHisCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomHisCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomHisCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomHisCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomHisCA, AtomHisCB, 110.6);
    addBondAngle(AtomHisCB, AtomHisCA, AtomC, 110.4);
    addBondAngle(AtomHisCA, AtomHisCB, AtomHisCG, 113.6);
    addBondAngle(AtomHisCB, AtomHisCG, AtomHisCD2, 130.8);
    addBondAngle(AtomHisCB, AtomHisCG, AtomHisND1, 123.2);
    addBondAngle(AtomHisCG, AtomHisCD2, AtomHisNE2, 109.2);
    addBondAngle(AtomHisCG, AtomHisND1, AtomHisCE1, 108.2);

    addBondLength(AtomC, AtomPheCA, 1.525);
    addBondLength(AtomNH1, AtomPheCA, 1.459);
    addBondLength(AtomPheCA, AtomPheCB, 1.535);
    addBondLength(AtomPheCB, AtomPheCG, 1.509);
    addBondLength(AtomPheCG, AtomPheCD1, 1.383);
    addBondLength(AtomPheCG, AtomPheCD2, 1.383);
    addBondLength(AtomPheCD1, AtomPheCE1, 1.388);
    addBondLength(AtomPheCD2, AtomPheCE2, 1.388);
    addBondLength(AtomPheCE2, AtomPheCZ, 1.369);

    addBondAngle(AtomNH1, AtomPheCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomPheCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomPheCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomPheCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomPheCA, AtomPheCB, 110.6); // NOT ENGH
    addBondAngle(AtomPheCB, AtomPheCA, AtomC, 110.1);
    addBondAngle(AtomPheCA, AtomPheCB, AtomPheCG, 113.9);
    addBondAngle(AtomPheCD1, AtomPheCG, AtomPheCD2, 118.299);
    addBondAngle(AtomPheCB, AtomPheCG, AtomPheCD2, 120.8);
    addBondAngle(AtomPheCB, AtomPheCG, AtomPheCD1, 120.8);
    addBondAngle(AtomPheCG, AtomPheCD2, AtomPheCE2, 120.8);
    addBondAngle(AtomPheCG, AtomPheCD1, AtomPheCE1, 120.8);
    addBondAngle(AtomPheCD2, AtomPheCE2, AtomPheCZ, 120.1);

    addBondLength(AtomC, AtomIleCA, 1.525);
    addBondLength(AtomNH1, AtomIleCA, 1.459);
    addBondLength(AtomIleCA, AtomIleCB, 1.544);
    addBondLength(AtomIleCB, AtomIleCG1, 1.536);
    addBondLength(AtomIleCG1, AtomIleCD1, 1.500);
    addBondLength(AtomIleCB, AtomIleCG2, 1.524);

    addBondAngle(AtomNH1, AtomIleCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomIleCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomIleCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomIleCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomIleCA, AtomIleCB, 110.8);
    addBondAngle(AtomIleCB, AtomIleCA, AtomC, 111.6);
    addBondAngle(AtomIleCA, AtomIleCB, AtomIleCG1, 110.0);
    addBondAngle(AtomIleCA, AtomIleCB, AtomIleCG2, 110.9);
    addBondAngle(AtomIleCG1, AtomIleCB, AtomIleCG2, 111.4);
    addBondAngle(AtomIleCB, AtomIleCG1, AtomIleCD1, 113.9);

    addBondLength(AtomC, AtomAlaCA, 1.525);
    addBondLength(AtomNH1, AtomAlaCA, 1.459);
    addBondLength(AtomAlaCA, AtomAlaCB, 1.520);

    addBondAngle(AtomNH1, AtomAlaCA, AtomC, C_CH1E_NH1_ANGLE);  // yes
    addBondAngle(AtomAlaCA, AtomC, AtomO, CH1E_C_O_ANGLE);  // yes
    addBondAngle(AtomC, AtomNH1, AtomAlaCA, C_NH1_CH1E_ANGLE); // no
    addBondAngle(AtomAlaCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomAlaCA, AtomAlaCB, 110.1);
    addBondAngle(AtomAlaCB, AtomAlaCA, AtomC, 110.1);

    addBondLength(AtomC, AtomLeuCA, 1.525);
    addBondLength(AtomNH1, AtomLeuCA, 1.459);
    addBondLength(AtomLeuCA, AtomLeuCB, 1.533);
    addBondLength(AtomLeuCB, AtomLeuCG, 1.521);
    addBondLength(AtomLeuCG, AtomLeuCD1, 1.524);
    addBondLength(AtomLeuCG, AtomLeuCD2, 1.524);

    addBondAngle(AtomNH1, AtomLeuCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomLeuCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomLeuCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomLeuCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomLeuCA, AtomLeuCB, 110.4);
    addBondAngle(AtomLeuCB, AtomLeuCA, AtomC, 110.2);
    addBondAngle(AtomLeuCA, AtomLeuCB, AtomLeuCG, 115.3);
    addBondAngle(AtomLeuCB, AtomLeuCG, AtomLeuCD1, 111.0);
    addBondAngle(AtomLeuCB, AtomLeuCG, AtomLeuCD2, 110.5);
    addBondAngle(AtomLeuCD1, AtomLeuCG, AtomLeuCD2, 110.5);

    addBondLength(AtomC, AtomTrpCA, 1.525);
    addBondLength(AtomNH1, AtomTrpCA, 1.459);
    addBondLength(AtomTrpCA, AtomTrpCB, 1.535);
    addBondLength(AtomTrpCB, AtomTrpCG, 1.498);
    addBondLength(AtomTrpCG, AtomTrpCD1, 1.363);
    addBondLength(AtomTrpCG, AtomTrpCD2, 1.432);
    addBondLength(AtomTrpCD1, AtomTrpNE1, 1.375);
    addBondLength(AtomTrpNE1, AtomTrpCE2, 1.371);
    addBondLength(AtomTrpCD2, AtomTrpCE2, 1.409);
    addBondLength(AtomTrpCD2, AtomTrpCE3, 1.399);
    addBondLength(AtomTrpCE2, AtomTrpCZ2, 1.393);
    addBondLength(AtomTrpCE3, AtomTrpCZ3, 1.382);
    addBondLength(AtomTrpCZ3, AtomTrpCH2, 1.369);
    addBondLength(AtomTrpCZ2, AtomTrpCH2, 1.396);

    addBondAngle(AtomNH1, AtomTrpCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomTrpCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomTrpCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomTrpCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomTrpCA, AtomTrpCB, 110.6);
    addBondAngle(AtomTrpCB, AtomTrpCA, AtomC, 110.4);
    addBondAngle(AtomTrpCA, AtomTrpCB, AtomTrpCG, 113.7);
    addBondAngle(AtomTrpCB, AtomTrpCG, AtomTrpCD1, 127.0);
    addBondAngle(AtomTrpCB, AtomTrpCG, AtomTrpCD2, 126.6);
    addBondAngle(AtomTrpCD1, AtomTrpCG, AtomTrpCD2, 106.3);
    addBondAngle(AtomTrpCG, AtomTrpCD1, AtomTrpNE1, 110.1);
    addBondAngle(AtomTrpCG, AtomTrpCD2, AtomTrpCE3, 133.9);
    addBondAngle(AtomTrpCD1, AtomTrpNE1, AtomTrpCE2, 109.0);
    addBondAngle(AtomTrpNE1, AtomTrpCE2, AtomTrpCZ2, 130.4);
    addBondAngle(AtomTrpCZ3, AtomTrpCE3, AtomTrpCD2, 118.8);
    addBondAngle(AtomTrpCD2, AtomTrpCE2, AtomTrpCZ2, 122.3);
    addBondAngle(AtomTrpCE3, AtomTrpCZ3, AtomTrpCH2, 121.2);

    addBondLength(AtomC, AtomTyrCA, 1.525);
    addBondLength(AtomNH1, AtomTyrCA, 1.459);
    addBondLength(AtomTyrCA, AtomTyrCB, 1.535);
    addBondLength(AtomTyrCB, AtomTyrCG, 1.512);
    addBondLength(AtomTyrCG, AtomTyrCD1, 1.387);
    addBondLength(AtomTyrCG, AtomTyrCD2, 1.387);
    addBondLength(AtomTyrCD1, AtomTyrCE1, 1.389);
    addBondLength(AtomTyrCD2, AtomTyrCE2, 1.389);
    addBondLength(AtomTyrCE2, AtomTyrCZ, 1.381);
    addBondLength(AtomTyrCZ, AtomTyrOH, 1.374);

    addBondAngle(AtomNH1, AtomTyrCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomTyrCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomTyrCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomTyrCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomTyrCA, AtomTyrCB, 110.5);
    addBondAngle(AtomC, AtomTyrCA, AtomTyrCB, 110.1);
    addBondAngle(AtomTyrCA, AtomTyrCB, AtomTyrCG, 113.9);
    addBondAngle(AtomTyrCB, AtomTyrCG, AtomTyrCD2, 121.0);
    addBondAngle(AtomTyrCD1, AtomTyrCG, AtomTyrCD2, 117.9);
    addBondAngle(AtomTyrCB, AtomTyrCG, AtomTyrCD1, 121.0);
    addBondAngle(AtomTyrCG, AtomTyrCD2, AtomTyrCE2, 121.3);
    addBondAngle(AtomTyrCG, AtomTyrCD1, AtomTyrCE1, 121.3);
    addBondAngle(AtomTyrCD2, AtomTyrCE2, AtomTyrCZ, 119.8);
    addBondAngle(AtomTyrCE2, AtomTyrCZ, AtomTyrOH, 120.1);

    addBondLength(AtomC, AtomAspCA, 1.525);
    addBondLength(AtomNH1, AtomAspCA, 1.459);
    addBondLength(AtomAspCA, AtomAspCB, 1.535);
    addBondLength(AtomAspCB, AtomAspCG, 1.513);
    addBondLength(AtomAspCG, AtomAspOD1, 1.249);
    addBondLength(AtomAspCG, AtomAspOD2, 1.249);

    addBondAngle(AtomNH1, AtomAspCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomAspCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomAspCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomAspCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomAspCA, AtomAspCB, 110.6);
    addBondAngle(AtomAspCB, AtomAspCA, AtomC, 110.4);
    addBondAngle(AtomAspCA, AtomAspCB, AtomAspCG, 113.4);
    addBondAngle(AtomAspCB, AtomAspCG, AtomAspOD1, 118.3);
    addBondAngle(AtomAspCB, AtomAspCG, AtomAspOD2, 118.3);
    addBondAngle(AtomAspOD1, AtomAspCG, AtomAspOD2, 123.4);

    addBondLength(AtomC, AtomThrCA, 1.525);
    addBondLength(AtomNH1, AtomThrCA, 1.459);
    addBondLength(AtomThrCA, AtomThrCB, 1.529);
    addBondLength(AtomThrCB, AtomThrOG1, 1.428);
    addBondLength(AtomThrCB, AtomThrCG2, 1.519);

    addBondAngle(AtomNH1, AtomThrCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomThrCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomThrCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomThrCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomThrCA, AtomThrCB, 110.3);
    addBondAngle(AtomThrCB, AtomThrCA, AtomC, 111.6);
    addBondAngle(AtomThrCA, AtomThrCB, AtomThrCG2, 112.4);
    addBondAngle(AtomThrCA, AtomThrCB, AtomThrOG1, 109.0);
    addBondAngle(AtomThrCG2, AtomThrCB, AtomThrOG1, 109.3);

    addBondLength(AtomC, AtomValCA, 1.525);
    addBondLength(AtomNH1, AtomValCA, 1.459);
    addBondLength(AtomValCA, AtomValCB, 1.543);
    addBondLength(AtomValCB, AtomValCG1, 1.524);
    addBondLength(AtomValCB, AtomValCG2, 1.524);

    addBondAngle(AtomNH1, AtomValCA, AtomC, C_CH1E_NH1_ANGLE);
    addBondAngle(AtomValCA, AtomC, AtomO, CH1E_C_O_ANGLE);
    addBondAngle(AtomC, AtomNH1, AtomValCA, C_NH1_CH1E_ANGLE);
    addBondAngle(AtomValCA, AtomC, AtomNH1, CH1E_C_NH1_ANGLE);
    addBondAngle(AtomNH1, AtomValCA, AtomValCB, 111.5);
    addBondAngle(AtomValCB, AtomValCA, AtomC, 111.4);
    addBondAngle(AtomValCA, AtomValCB, AtomValCG1, 110.9);
    addBondAngle(AtomValCA, AtomValCB, AtomValCG2, 110.9);
    addBondAngle(AtomValCG1, AtomValCB, AtomValCG2, 110.9);

}

std::string GeomTable::getResCode(std::string three)
{
    GeomTable table = getGeomTable();

    if (table._three2OneCode.count(three))
    {
        return table._three2OneCode[three];
    }

    return "?";
}
