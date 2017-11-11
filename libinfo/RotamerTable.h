//
//  RotamerTable.hpp
//  vagabond
//
//  Created by Helen Ginn on 10/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef RotamerTable_hpp
#define RotamerTable_hpp

#include <string>
#include <stdio.h>
#include <vector>
#include <map>

typedef struct
{
	std::string firstAtom;
	std::string secondAtom;
	double torsion;
} TorsionAngle;

typedef struct
{
	std::string name;
	std::vector<TorsionAngle> torsions;
	double allOccupancy;
	double alphaOccupancy;
	double betaOccupancy;
	double otherOccupancy;
} Rotamer;

typedef std::map<std::string, std::vector<Rotamer> > ResidueRotamerMap;

class RotamerTable
{
public:
	RotamerTable();

	std::vector<Rotamer> rotamersFor(std::string residue);

	static RotamerTable *getTable()
	{
		return &_rotamerTable;
	}
	
private:
	static RotamerTable _rotamerTable;

	ResidueRotamerMap _residueRotamerMap;

	static TorsionAngle createTorsionAngle(std::string atom1, std::string atom2,
										   double value);
	static void setOccupancies(Rotamer *rot, double all, double alpha,
							   double beta, double other);

	void addRotamerToTable(std::string resName, Rotamer rot);
};

#endif /* RotamerTable_hpp */
