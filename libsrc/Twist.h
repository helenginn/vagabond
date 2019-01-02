// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
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

#ifndef __vagabond__Twist_h
#define __vagabond__Twist_h

#include "shared_ptrs.h"
#include "Parser.h"
#include "ExplicitModel.h"
#include "Whack.h"

/**
 * \class Twist
 * \brief Provides an anchor-independent twist around a bond which produces
 * the ability to change the torsion angle backwards.
 */

class Twist : public boost::enable_shared_from_this<Twist>
{
public:
	Twist();

	/** Add the Whack to an anchor point */
	void addToAppliedModel(ExplicitModelPtr applied);

	/** Bond to which twist is applied upstream/downstream */
	void setBond(BondPtr bond);
	
	/** To be used by AnchorPtr */
	void disable()
	{
		_enabled = false;
	}
	
	/** To be used by AnchorPtr */
	void enable()
	{
		_enabled = true;
	}

	/** To be called by an Anchor object to modify its sample positions. */
	void applyToAnchorSamples(std::vector<BondSample> &anchSamp);

	std::string getClassName()
	{
		return "Twist";
	}

	/* Returns a pointer to the anchor */
	ExplicitModelPtr getAppliedModel()
	{
		return _applied.lock();
	}
	
	/* Returns twist angle in radians */
	double getTorsionCorrection()
	{
		if (!_enabled) return 0;
		return -_twist;
	}
	
	static double getTwist(void *object)
	{
		return static_cast<Twist *>(object)->_twist;
	}
	
	static void setTwist(void *object, double val);

	void saveSamples();
private:
	bool _enabled;
	bool _valid;
	double _twist;
	
	std::vector<BondSample> _samples;
	ExplicitModelWkr _applied;
	BondPtr _bond;
};

#endif
