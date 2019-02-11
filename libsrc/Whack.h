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
#include "Parser.h"

/**
 * \class Whack
 * \brief Whack classes allow for a bi-directional kick to a bond to
 * escape dependency on the anchor point.
 *
 * Whacks allow for a kick to be applied to a bond with a proportion
 * going forwards and the rest going backwards. It's a bit stronger than
 * a kick - so it's called a Whack.
 */

inline vec3 rotate_round_bond(vec3 start, vec3 centre, mat3x3 rot)
{
	vec3_subtract_from_vec3(&start, centre);
	mat3x3_mult_vec(rot, &start);
	vec3_add_to_vec3(&start, centre);

	return start;
}


class Whack : public Parser
{
public:
	Whack();
	
	WhackPtr shared_from_this()
	{
		return ToWhackPtr(Parser::shared_from_this());
	}
	
	/** Add the Whack to an anchor point to make it valid, otherwise
	 * will shout at Helen. */
	void addToAnchor(AnchorPtr anchor);
	
	/** Bond to which kick is applied upstream/downstream */
	void setBond(BondPtr bond);
	
	/** Return bond to which Whack is attached */
	BondPtr getBond()
	{
		return _bond;
	}
	
	/** Set the proportion of the kick which should be passed into the
	 * downstream bonds. 1 = like a normal kick. 0 = backwards, towards
	 * the anchor point, only. */
	static void setKick(void *object, double kick);

	/** Set the magnitude of the kick which is to be shared between
	 * the forwards and backwards directions */
	static void setWhack(void *object, double whack);
	
	/** Get the magnitude of the whack */
	static double getWhack(void *object)
	{
		return static_cast<Whack *>(object)->_whack;
	}

	/** Get the magnitude of the kick */
	static double getKick(void *object)
	{
		return static_cast<Whack *>(object)->_kick;
	}
	
	/** To be called by an Anchor object to modify its sample positions. */
	void applyToAnchorSamples(std::vector<BondSample> &anchSamp);

	/** Refresh information used to calculate the Whack. */
	void applyKick();
	
	/** Save samples from _bond. */
	void saveSamples();
	
	AnchorPtr getAnchor()
	{
		return _anchor.lock();
	}
	
	/* Whack is invalid if it is attached to a bond for which the child is
	 * incapable of refining flexibility. */
	bool isValid()
	{
		return _valid;
	}
	
	void disable()
	{
		_enabled = false;
		applyKick();
	}
	
	void enable()
	{
		_enabled = true;
		applyKick();
	}
	
	virtual std::string getClassName()
	{
		return "Whack";
	}

	/** If the number of samples in the ensemble has changed, returns true */
	bool needsRefresh(std::vector<BondSample> &anchSamp);
protected:
	virtual std::string getParserIdentifier();
	
	virtual void addProperties();
	virtual void linkReference(ParserPtr object, std::string category);
	virtual void postParseTidy();
	
private:
	std::vector<BondSample> _samples;
	double _kick;
	double _whack;
	bool _valid;
	bool _enabled;
	
	AnchorWkr _anchor;
	BondPtr _bond;
};


#endif
