// Vagabond
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#ifndef __Vagabond__Angler__
#define __Vagabond__Angler__

#include "shared_ptrs.h"
#include "Parser.h"

typedef enum
{
	AngleCNA,
	AngleNAB,
	AngleNAC,
	AngleBAC,
	AngleACN,
	AngleACO,
	AngleOCN,
	AngleNone
} AngleType;

typedef enum
{
	ResiNonPGIV,
	ResiIleVal,
	ResiProline,
	ResiGlycine,
} ResiType;

class Angler : public Parser
{
public:
	AnglerPtr shared_from_this()
	{
		return ToAnglerPtr(Parser::shared_from_this());
	}

	Angler();
	
	virtual ~Angler() {};
	
	/* returns true if valid */
	bool setupTable();
	
	double getMainAngle();
	double getSisterAngle();
	double getOtherAngle();
	
	void setNextIsProline(bool pro)
	{
		_nextIsPro = pro;
	}
	
	void setBonds(BondPtr angle, BondPtr phi, BondPtr psi)
	{
		_angledBond = angle;
		_phiBond = phi;
		_psiBond = psi;
	}
	
	void setAngledBond(BondPtr bond)
	{
		_angledBond = bond;
	}

	void setPhiTorsionBond(BondPtr bond)
	{
		_phiBond = bond;
	}

	void setPsiTorsionBond(BondPtr bond)
	{
		_psiBond = bond;
	}

	virtual std::string getClassName()
	{
		return "Angler";
	}
protected:
	virtual std::string getParserIdentifier();
	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void postParseTidy();

private:
	AngleType mainAngleType(AtomPtr earlier, AtomPtr major,
	                        AtomPtr sister, AtomPtr minor);
	void assignTable(ResiType res, AngleType ang, double **where);
	double getAngle(bool report = false, double *which = NULL);
	
	bool _nextIsPro;
	BondPtr _angledBond;
	BondPtr _phiBond;
	BondPtr _psiBond;

	double *_table;
	double *_sisterTable;
	double *_otherTable;
};

#endif
