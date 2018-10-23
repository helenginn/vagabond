// Vagabond : bond-based macromolecular model refinement
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

#ifndef __vagabond__Whack_h__
#define __vagabond__Whack_h__

#include "shared_ptrs.h"
#include "ExplicitModel.h"

/**
 * \class Whack
 * \brief Whack classes allow for a bi-directional kick to a bond to
 * escape dependency on the anchor point.
 *
 * Whacks allow for a kick to be applied to a bond with a proportion
 * going forwards and the rest going backwards. It's a bit stronger than
 * a kick - so it's called a Whack.
 */

class Whack : public boost::enable_shared_from_this<Whack>
{
public:
	Whack();
	
	/** Add the Whack to an anchor point to make it valid, otherwise
	 * will shout at Helen. */
	void addToAnchor(AnchorPtr anchor);
	
	/** Bond to which kick is applied upstream/downstream */
	void setBond(BondPtr bond);
	
	/** Set the proportion of the kick which should be passed into the
	 * downstream bonds. 1 = like a normal kick. 0 = backwards, towards
	 * the anchor point, only. */
	static void setKick(void *object, double kick);

	/** Set the magnitude of the kick which is to be shared between
	 * the forwards and backwards directions */
	static void setWhack(void *object, double whack);
	
	/** To be called by an Anchor object to modify its sample positions. */
	void applyToAnchorSamples(std::vector<BondSample> &anchSamp);

	/** Refresh information used to calculate the Whack. */
	void applyKick();
private:
	std::vector<BondSample> _samples;
	double _kick;
	double _whack;
	
	AnchorWkr _anchor;
	BondPtr _bond;
};


#endif
