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

#ifndef __vagabond__Motion__
#define __vagabond__Motion__

#include "RefineMat3x3.h"
#include "ExplicitModel.h"

class FlexGlobal;

/**
 * \class Motion
 * \brief Whole-molecule motions which can be applied to multiple anchor
 * points if needed. 
 **/

class Motion : public Parser
{
public:
	Motion();
	~Motion();
	void addToPolymer(PolymerPtr pol);
	void removeFromPolymer(PolymerPtr pol);
	void updateAtoms();

	void applyTranslations(std::vector<BondSample> &stored,
	                       bool isomorphous = false);
	void applyRotations(std::vector<BondSample> &stored);
	void applyMotions(std::vector<BondSample> &stored);

	void addTranslationParameters(RefinementStrategyPtr strategy);
	void refine(bool reciprocal = false);
	void rigidRefine();
	void removeAtom(AtomPtr atom);
	
	int librationCount()
	{
		return _quats.size();
	}
	
	size_t screwCount()
	{
		return _screws.size();
	}

	void addLibrationParameters(RefinementStrategyPtr strategy,
	                              int num = 0);
	void addScrewParameters(RefinementStrategyPtr strategy,
	                              int num = -1);
	void deleteLastScrew();
	
	void reset();
	static void setScale(void *object, double scale);
	static double getScale(void *object);

	static void setTransScale(void *object, double scale)
	{
		static_cast<Motion *>(object)->_transScale = scale;
	}

	static double getTransScale(void *object)
	{
		return static_cast<Motion *>(object)->_transScale;

	}
	
	void absorbScale();

	MotionPtr shared_from_this()
	{
		return ToMotionPtr(BaseParser::shared_from_this());
	}
	
	virtual std::string getClassName()
	{
		return "Motion";
	}
	
	virtual std::string getParserIdentifier()
	{
		return "Motion_" + _name;
	}
	
	void setName(std::string name)
	{
		_name = name;
	}
	
	std::string getName()
	{
		return _name;
	}

	bool hasRefined()
	{
		return _refined;
	}

	void attachTargetToRefinement(RefinementStrategyPtr strategy,
	                              FlexGlobal &target, bool recip = false);
protected:
	virtual void addProperties();
	virtual void postParseTidy();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void linkReference(BaseParserPtr object, std::string category);
private:
	void deleteQuats();

	mat3x3 getOverallRotation();

	RefineMat3x3Ptr _trans;
	std::string _name;

	/* tmp for loading in/out of a vbond file */
	std::vector<vec3> _tmpQuats;
	std::vector<vec3> _tmpScrews;

	/* contains rotation information */
	std::vector<Quat4Refine *> _quats;
	
	std::vector<Quat4Refine *> _screws;

	Quat4Refine *_rotation;
	vec3 _r;
	Quat4Refine *_displacement;
	vec3 _d;

	double _scale;
	double _transScale;
	bool _refined;
	vec3 _centre;
	AtomGroupPtr _allAtoms;
	AtomGroupPtr _allBackbone;
	
	std::vector<MoleculeWkr> _molecules;
};

#endif
