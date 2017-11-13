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


}
