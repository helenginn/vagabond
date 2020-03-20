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

#ifndef __vagabond__Sponge__
#define __vagabond__Sponge__

#include "Novalent.h"
#include "MapScoreWorkspace.h"

typedef std::map<AtomWkr, double> Restraint;

class Sponge : public Novalent
{
public:
	SpongePtr shared_from_this()
	{
		return ToSpongePtr(Parser::shared_from_this());
	}

	Sponge(AtomPtr water);
	Sponge();
	virtual ~Sponge();

	virtual std::vector<BondSample> *getManyPositions(void *object = NULL);

	std::string getClassName()
	{
		return "Sponge";
	}
	
	AtomPtr getConnectedAtom(int i);
	
	virtual bool isDisabled()
	{
		return _disabled;
	}
	
	double distance(int i)
	{
		return _restraints[getConnectedAtom(i)];
	}
	
	size_t connectedCount();
	
	bool hasConnectedAtom(AtomPtr atom);
	
	void setActive(int n)
	{
		_n = n;
	}
	
	virtual std::string shortDesc();
	
	void singleRefine(bool others = false);
	void setupRestraints();
	void findConnections();
	void initialConnections();
	void copyActiveToFinalPos();
	void addRestraintsToStrategy(RefinementStrategyPtr strategy);

	void recalculated()
	{
		_changedSamples = false;
		_recalcFinal = false;
	}
	
	bool needsRecalc()
	{
		return (_changedSamples || _recalcFinal);
	}
	
	static void setX(void *object, double value)
	{
		Sponge *s = static_cast<Sponge *>(object);
		s->_storedSamples[s->_n].start.x = value;
	}
	
	static void setY(void *object, double value)
	{
		Sponge *s = static_cast<Sponge *>(object);
		s->_storedSamples[s->_n].start.y = value;
	}
	
	static void setZ(void *object, double value)
	{
		Sponge *s = static_cast<Sponge *>(object);
		s->_storedSamples[s->_n].start.z = value;
	}
	
	static double getX(void *object)
	{
		Sponge *s = static_cast<Sponge *>(object);
		return s->_storedSamples[s->_n].start.x;
	}
	
	static double getY(void *object)
	{
		Sponge *s = static_cast<Sponge *>(object);
		return s->_storedSamples[s->_n].start.y;
	}
	
	static double getZ(void *object)
	{
		Sponge *s = static_cast<Sponge *>(object);
		return s->_storedSamples[s->_n].start.z;
	}
private:
	void randomConnections();
	AtomGroupPtr correlGroup();
	WaterNetworkPtr getNetwork();

	AtomGroupPtr _close;
	MapScoreWorkspace _workspace;
	std::vector<AtomPtr> _candidates;

	double _preScore;
	bool _disabled;
	int _n;
	Restraint _restraints;
};

#endif
