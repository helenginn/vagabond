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

typedef enum
{
	AtomUnassigned,
	AtomC,
	AtomCH1E,
	AtomCH2E,
	AtomCH3E,
	AtomNH3,
	AtomOH1,
} AtomType;

typedef std::pair<AtomType, AtomType> AtomPair;
typedef std::pair<AtomPair, AtomType> AtomTrio;

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

private:
	void addBondAngle(AtomType atom1, AtomType atom2,
					  AtomType atom3, double angle);

	void addBondLength(AtomType atom1, AtomType atom2, double length);
	static GeomTable _geomTable;

	std::map<AtomPair, double> _bondLengths;
	std::map<AtomTrio, double> _bondAngles;
};



#endif /* defined(__vagabond__GeomTable__) */
