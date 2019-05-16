//
//  GeomTable.h
//  vagabond
//
//  Created by Helen Ginn on 07/08/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__GeomTable__
#define __vagabond__GeomTable__

#include <stdio.h>
#include <map>
#include <string>

typedef enum
{
	AtomUnassigned,
        AtomUnknown,
	AtomC,
	AtomCH1E,
	AtomCH2E,
	AtomCH3E,
	AtomC5,
	AtomCRH,
	AtomCR1H,
	AtomCR1E,
	AtomNH1,
	AtomNH3,
	AtomOH1,
	AtomO,
	AtomOXT,
	AtomSM,

	AtomGlyCA,
	AtomAlaCA, AtomAlaCB,
	AtomLeuCA, AtomLeuCB, AtomLeuCG, AtomLeuCD1, AtomLeuCD2,
	AtomPheCA, AtomPheCB, AtomPheCG, AtomPheCD1, AtomPheCD2, AtomPheCE1, AtomPheCE2, AtomPheCZ,
	AtomTyrCA, AtomTyrCB, AtomTyrCG, AtomTyrCD1, AtomTyrCD2, AtomTyrCE1, AtomTyrCE2, AtomTyrCZ, AtomTyrOH,
	AtomHisCA, AtomHisCB, AtomHisCG, AtomHisCE1, AtomHisCD2, AtomHisND1, AtomHisNE2,
	AtomMetCA, AtomMetCB, AtomMetCG, AtomMetSD, AtomMetCE,
	AtomLysCA, AtomLysCB, AtomLysCG, AtomLysCD, AtomLysCE, AtomLysNZ,
	AtomSerCA, AtomSerCB, AtomSerOG,
	AtomCysCA, AtomCysCB, AtomCysCBS, AtomCysSG, AtomCysSGS,
	AtomThrCA, AtomThrCB, AtomThrOG1, AtomThrCG2,
	AtomValCA, AtomValCB, AtomValCG1, AtomValCG2,
	AtomAspCA, AtomAspCB, AtomAspCG, AtomAspOD1, AtomAspOD2,
	AtomAsnCA, AtomAsnCB, AtomAsnCG, AtomAsnOD1, AtomAsnND2,
	AtomGlnCA, AtomGlnCB, AtomGlnCG, AtomGlnCD, AtomGlnOE1, AtomGlnNE2,
	AtomGluCA, AtomGluCB, AtomGluCG, AtomGluCD, AtomGluOE1, AtomGluOE2,
	AtomIleCA, AtomIleCB, AtomIleCG1, AtomIleCG2, AtomIleCD1,
	AtomProNH1, AtomProC, AtomProCA, AtomProCB, AtomProCG, AtomProCD,
	AtomArgCA, AtomArgCB, AtomArgCG, AtomArgCD, AtomArgNE, AtomArgCZ, AtomArgNH1, AtomArgNH2,
	AtomTrpCA, AtomTrpCB, AtomTrpCG, AtomTrpCD1, AtomTrpCD2, AtomTrpNE1, AtomTrpCE2, AtomTrpCE3, AtomTrpCZ2, AtomTrpCZ3, AtomTrpCH2,

} AtomType;

typedef std::pair<AtomType, AtomType> AtomPair;
typedef std::pair<AtomPair, AtomType> AtomTrio;
typedef std::pair<std::string, std::string> AtomIdentity;

class GeomTable
{
public:
	static GeomTable getGeomTable()
	{
		return _geomTable;
	}

	GeomTable();
	double getBondLength(AtomType atom1, AtomType atom2);
	double getBondAngle(AtomType atom1, AtomType atom2, AtomType atom3);
	int getChiralCentre(AtomType atom1, AtomType atom2, 
                                  AtomType atom3);
	AtomType getType(std::string, std::string);

	static std::string getResCode(std::string threeLetter);

private:
	void addBondAngle(AtomType atom1, AtomType atom2,
					  AtomType atom3, double angle);
	void addIdentityToType(std::string, std::string, AtomType type);
	void addBondLength(AtomType atom1, AtomType atom2, double length);
	void addChiralCentre(AtomType atom1, AtomType atom2, 
	                     AtomType atom3, int sign);
	void addSingleChiral(AtomType atom1, AtomType atom2, 
	                     AtomType atom3, int sign);

	static GeomTable _geomTable;
	std::map<std::string, std::string> _three2OneCode;

	std::map<AtomPair, double> _bondLengths;
	std::map<AtomTrio, double> _bondAngles;
	std::map<AtomTrio, int> _chirals;
	std::map<AtomIdentity, AtomType> _identityToType;
};



#endif /* defined(__vagabond__GeomTable__) */
