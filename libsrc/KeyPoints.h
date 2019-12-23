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

#ifndef __Vagabond__KeyPoints__h
#define __Vagabond__KeyPoints__h

#include <vector>
#include "shared_ptrs.h"
#include "FlexGlobal.h"
#include "Param.h"
#include "Parser.h"

typedef struct
{
	Param res;
	Param phi;
	Param psi;
	Param kick;
} WayPoint;

typedef enum
{
	WayPointPhi,
	WayPointPsi,
	WayPointKick
} WayPointType;

/** 
 * \class KeyPoints
 * \brief specify key points for a given monomer to derive phi/psi
 * modifications for a Bond */

class KeyPoints : public Parser
{
public:
	KeyPointsPtr shared_from_this()
	{
		return ToKeyGroupPtr(Parser::shared_from_this());
	}

	KeyPoints();
	
	void setPolymer(PolymerPtr polymer);

	double getContribution(BondPtr bond, WayPointType type);
	double getPhiContribution(BondPtr bond);
	double getPsiContribution(BondPtr bond);
	double getKickContribution(BondPtr bond);
	
	bool refineKeyPoints();
protected:
	virtual std::string getClassName()
	{
		return "KeyPoints";
	}

	/* should only be one identifier per polymer */
	virtual std::string getParserIdentifier()
	{
		return "KeyPoints";
	}
	
	
	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void linkReference(BaseParserPtr object, std::string category);
	virtual void postParseTidy();
private:
	static double score(void *object);
	std::vector<WayPoint> _points;

	PolymerPtr _polymer;
	FlexGlobal _global;
	std::vector<vec3> _tmpWays;
	std::vector<vec3> _tmpKicks;
};

#endif
