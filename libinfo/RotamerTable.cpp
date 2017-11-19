//
//  RotamerTable.cpp
//  vagabond
//
//  Created by Helen Ginn on 10/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "RotamerTable.h"
#include <algorithm>

RotamerTable RotamerTable::_rotamerTable = RotamerTable();

TorsionAngle RotamerTable::createTorsionAngle(std::string atom1, std::string atom2,
											  double value)
{
	TorsionAngle angle;
	angle.firstAtom = atom1;
	angle.secondAtom = atom2;
	angle.torsion = value;

	return angle;
}

void RotamerTable::setOccupancies(Rotamer *rot, double all, double alpha,
								  double beta, double other)
{
	rot->allOccupancy = all;
	rot->alphaOccupancy = alpha;
	rot->betaOccupancy = beta;
	rot->otherOccupancy = other;
}

void RotamerTable::addRotamerToTable(std::string resName, Rotamer rot)
{
	_residueRotamerMap[resName].push_back(rot);
}

bool rotamerMoreCommon(Rotamer rot1, Rotamer rot2)
{
	return (rot1.allOccupancy > rot2.allOccupancy);
}

std::vector<Rotamer> RotamerTable::rotamersFor(std::string residue)
{
	if (!_residueRotamerMap.count(residue))
	{
		return std::vector<Rotamer>();
	}

	std::vector<Rotamer> rotamers = _residueRotamerMap[residue];
	std::sort(rotamers.begin(), rotamers.end(), rotamerMoreCommon);

	return rotamers;
}

RotamerTable::RotamerTable()
{
	/** ASPARAGINE / ASN / N **/
	{
		Rotamer rotamer;
		rotamer.name = "p-10";
		setOccupancies(&rotamer, 0.07, 0.00, 0.01, 0.10);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", 62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -10.0));
		addRotamerToTable("asn", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "p30";
		setOccupancies(&rotamer, 0.09, 0.005, 0.07, 0.12);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", 62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 30.0));
		addRotamerToTable("asn", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t-20";
		setOccupancies(&rotamer, 0.12, 0.05, 0.21, 0.12);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -174.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -20.0));
		addRotamerToTable("asn", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t30";
		setOccupancies(&rotamer, 0.15, 0.13, 0.18, 0.15);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 30.0));
		addRotamerToTable("asn", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m-20";
		setOccupancies(&rotamer, 0.39, 0.65, 0.28, 0.33);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -20.0));
		addRotamerToTable("asn", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m-80";
		setOccupancies(&rotamer, 0.08, 0.08, 0.09, 0.08);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -75.0));
		addRotamerToTable("asn", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m120";
		setOccupancies(&rotamer, 0.04, 0.03, 0.03, 0.04);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 120.0));
		addRotamerToTable("asn", rotamer);
	}

	/** ASPARTATE / ASP / D **/

	{
		Rotamer rotamer;
		rotamer.name = "p-10";
		setOccupancies(&rotamer, 0.10, 0.01, 0.02, 0.13);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", 62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -10.0));
		addRotamerToTable("asp", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "p30";
		setOccupancies(&rotamer, 0.09, 0.01, 0.05, 0.12);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", 62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 30.0));
		addRotamerToTable("asp", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t0";
		setOccupancies(&rotamer, 0.21, 0.08, 0.44, 0.20);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 0.0));
		addRotamerToTable("asp", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t70";
		setOccupancies(&rotamer, 0.06, 0.11, 0.07, 0.04);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 65.0));
		addRotamerToTable("asp", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m-20";
		setOccupancies(&rotamer, 0.51, 0.77, 0.38, 0.47);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -70.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -15.0));
		addRotamerToTable("asp", rotamer);
	}

	/** LYSINE / LYS / K **/

	{
		Rotamer rotamer;
		rotamer.name = "ptpt";
		setOccupancies(&rotamer, 0.01, 0.00, 0.02, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE", 180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "pttp";
		setOccupancies(&rotamer, 0.01, 0.00, 0.01, 0.02);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "pttt";
		setOccupancies(&rotamer, 0.02, 0.00, 0.04, 0.03);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE", 180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "pttm";
		setOccupancies(&rotamer, 0.01, 0.00, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE", -65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "ptmt";
		setOccupancies(&rotamer, 0.01, 0.00, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", 180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD", -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE", 180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tptp";
		setOccupancies(&rotamer, 0.01, 0.01, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",   65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tptt";
		setOccupancies(&rotamer, 0.03, 0.05, 0.01, 0.02);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tptm";
		setOccupancies(&rotamer, 0.01, 0.01, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  -65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "ttpp";
		setOccupancies(&rotamer, 0.01, 0.01, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",   65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "ttpt";
		setOccupancies(&rotamer, 0.02, 0.02, 0.05, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tttp";
		setOccupancies(&rotamer, 0.04, 0.05, 0.05, 0.03);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",   65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tttt";
		setOccupancies(&rotamer, 0.13, 0.17, 0.19, 0.10);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tttm";
		setOccupancies(&rotamer, 0.03, 0.04, 0.02, 0.03);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  -65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "ttmt";
		setOccupancies(&rotamer, 0.02, 0.02, 0.04, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "ttmm";
		setOccupancies(&rotamer, 0.01, 0.01, 0.00, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  -65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mptt";
		setOccupancies(&rotamer, 0.01, 0.00, 0.00, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -90.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mtpp";
		setOccupancies(&rotamer, 0.01, 0.01, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",   65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mtpt";
		setOccupancies(&rotamer, 0.03, 0.04, 0.02, 0.03);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",   68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mttp";
		setOccupancies(&rotamer, 0.03, 0.02, 0.04, 0.04);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",   65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mttt";
		setOccupancies(&rotamer, 0.20, 0.23, 0.14, 0.21);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mttm";
		setOccupancies(&rotamer, 0.05, 0.03, 0.05, 0.06);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  -65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mtmt";
		setOccupancies(&rotamer, 0.03, 0.06, 0.02, 0.03);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mtmm";
		setOccupancies(&rotamer, 0.01, 0.00, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  -65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mmtp";
		setOccupancies(&rotamer, 0.01, 0.01, 0.00, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",   65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mmtt";
		setOccupancies(&rotamer, 0.06, 0.03, 0.05, 0.08);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mmtm";
		setOccupancies(&rotamer, 0.01, 0.01, 0.01, 0.02);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  -65.0));
		addRotamerToTable("lys", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mmmt";
		setOccupancies(&rotamer, 0.01, 0.01, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -68.0));
		rotamer.torsions.push_back(createTorsionAngle("CD", "CE",  180.0));
		addRotamerToTable("lys", rotamer);
	}

	/* SERINE / SER / S */

	{
		Rotamer rotamer;
		rotamer.name = "p";
		setOccupancies(&rotamer, 0.48, 0.33, 0.36, 0.55);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		addRotamerToTable("ser", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t";
		setOccupancies(&rotamer, 0.22, 0.22, 0.34, 0.18);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		addRotamerToTable("ser", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m";
		setOccupancies(&rotamer, 0.29, 0.44, 0.29, 0.25);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -65.0));
		addRotamerToTable("ser", rotamer);
	}

	/* HISTIDINE / HIS / H */

	{
		Rotamer rotamer;
		rotamer.name = "p-80";
		setOccupancies(&rotamer, 0.09, 0.00, 0.06, 0.13);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -75.0));
		addRotamerToTable("his", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "p80";
		setOccupancies(&rotamer, 0.04, 0.00, 0.04, 0.06);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  80.0));
		addRotamerToTable("his", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t-160";
		setOccupancies(&rotamer, 0.05, 0.05, 0.14, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG", -165.0));
		addRotamerToTable("his", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t-80";
		setOccupancies(&rotamer, 0.11, 0.17, 0.09, 0.09);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -80.0));
		addRotamerToTable("his", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "t60";
		setOccupancies(&rotamer, 0.16, 0.24, 0.17, 0.12);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   60.0));
		addRotamerToTable("his", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m-70";
		setOccupancies(&rotamer, 0.29, 0.26, 0.30, 0.30);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -70.0));
		addRotamerToTable("his", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m170";
		setOccupancies(&rotamer, 0.07, 0.09, 0.03, 0.09);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  165.0));
		addRotamerToTable("his", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "m80";
		setOccupancies(&rotamer, 0.13, 0.14, 0.10, 0.14);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   85.0));
		addRotamerToTable("his", rotamer);
	}

	/* GLUTAMATE / GLU / E */

	{
		Rotamer rotamer;
		rotamer.name = "pt-20";
		setOccupancies(&rotamer, 0.05, 0.01, 0.09, 0.07);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",   62.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -20.0));
		addRotamerToTable("glu", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "pm0";
		setOccupancies(&rotamer, 0.02, 0.00, 0.00, 0.04);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",   70.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -80.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",    0.0));
		addRotamerToTable("glu", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tp10";
		setOccupancies(&rotamer, 0.06, 0.10, 0.02, 0.06);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   65.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",   10.0));
		addRotamerToTable("glu", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tt0";
		setOccupancies(&rotamer, 0.24, 0.25, 0.42, 0.18);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",   00.0));
		addRotamerToTable("glu", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "tm-20";
		setOccupancies(&rotamer, 0.01, 0.01, 0.01, 0.01);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB", -177.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -80.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -25.0));
		addRotamerToTable("glu", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mp0";
		setOccupancies(&rotamer, 0.06, 0.01, 0.02, 0.10);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",   85.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",    0.0));
		addRotamerToTable("glu", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mt-10";
		setOccupancies(&rotamer, 0.33, 0.36, 0.29, 0.32);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -67.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  180.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -10.0));
		addRotamerToTable("glu", rotamer);
	}

	{
		Rotamer rotamer;
		rotamer.name = "mm-40";
		setOccupancies(&rotamer, 0.13, 0.19, 0.07, 0.12);
		rotamer.torsions.push_back(createTorsionAngle("CA", "CB",  -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CB", "CG",  -65.0));
		rotamer.torsions.push_back(createTorsionAngle("CG", "CD",  -40.0));
		addRotamerToTable("glu", rotamer);
	}
}
