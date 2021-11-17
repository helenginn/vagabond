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

#ifndef __vagabond__Polymer__
#define __vagabond__Polymer__

#include <stdio.h>
#include "shared_ptrs.h"
#include "Molecule.h"
#include <vector>
#include <map>
#include "Options.h"

/**
 * \class Polymer
 * \brief A subclass of Molecule which contains a series of Monomer objects,
 * forming a polymer chain.
 */

class Polymer : public Molecule
{
public:
	Polymer();
	virtual ~Polymer() {}

	void closenessSummary();
	void addMonomer(MonomerPtr monomer);
	virtual void summary();
	virtual void tieAtomsUp();
	virtual void refine(CrystalPtr target, RefinementType rType);
	virtual void setStream(std::ostream *str);
	
	static void refinePositions(CrystalPtr cryst, PolymerPtr pol, int total);
	
	virtual void graph(std::string graphName);

	void scaleSidechainsToBFactor();
	void refineBackbone();
	
	virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal,
	                            int conformer = -1);

	virtual void reportParameters();
	void downWeightResidues(int start, int end, double value);

	void applyPolymerChanges();
	void refineToEnd(int monNum, CrystalPtr target, RefinementType rType);
	double refineRange(int start, int end, 
	                 CrystalPtr target, RefinementType rType);

	bool test();
	AnchorPtr getAnchorModel();
	void findAnchorNearestCentroid();
	void hydrogenateContents();
	void checkChainContinuity();

	bool isWhacking();
	int _whacked;
	void whack();
	void whackMonomer(MonomerPtr mon);
	AtomGroupPtr monomerRange(int start, int end, bool side = false);
	void ramachandranPlot();
	virtual void removeAtom(AtomPtr atom);

	void setAnchor(int num)
	{
		_anchorNum = num;
	}

	int getAnchor()
	{
		return _anchorNum;
	}

	MonomerPtr getMonomer(int i)
	{
		if (_monomers.count(i))
		{
			return _monomers[i];
		}

		return MonomerPtr();
	}

	int monomerBegin()
	{
		if (_monomers.size() == 0)
		{
			return 0;
		}

		return _monomers.begin()->first;
	}

	int monomerEnd()
	{
		if (_monomers.size() == 0)
		{
			return 0;
		}

		std::map<long, MonomerPtr>::iterator it = _monomers.end();
		it--;
		return it->first + 1;
	}
	
	int monomerCount()
	{
		//        return _monomers.size();
		return _totalMonomers;
	}

	bool hasResidue(int resNum);

	virtual std::string getClassName()
	{
		return "Polymer";
	}

	PolymerPtr shared_from_this()
	{
		return ToPolymerPtr(Molecule::shared_from_this());
	}

	void resetMotion();
	void refineMotions();
	void redefineMotion();
	void resetSidechains();
	void removeIntramolecularMotion();
	bool refineLocalFlexibility(bool magic = false);
	virtual void addProperties();
	virtual void addObject(ParserPtr object, std::string category);
	virtual void postParseTidy();
	
	AtomGroupPtr getAllBackbone();
	
	std::map<long, MonomerPtr>::iterator beginMonomer()
	{
		return _monomers.begin();
	}
	
	double getKickShift()
	{
		return _kickShift;
	}
	
	void setKickShift(double shift)
	{
		_kickShift = shift;
	}
	
	std::string getGraphName()
	{
		return _graphName;
	}
	
	virtual void addParamCounts(int *pos, int *flex)
	{
		*pos += _positionalParams;
		*flex += _flexibilityParams;
	}

	void scoreMonomers();
	
	void rigAnchor();
protected:
	virtual double getScore()
	{
		propagateChange();
		return Sampler::getScore();
	}

private:
	void setupKeyPoints();
	void refineMonomer(MonomerPtr monomer, CrystalPtr target,
	                   RefinementType rType);

	std::map<long, MonomerPtr> _monomers;

	int _anchorNum;
	double _startB;
	double _kick;
	double _kickShift;
	int _totalMonomers;
	int _flexibilityParams;
	int _positionalParams;
	std::string _graphName;

	std::map<MonomerPtr, vec3> _scores;
	double _fullScore;
	KeyPointsPtr _keyPoints;

	AnchorPtr _anchor;
	AtomGroupPtr _allBackbones;

};

#endif /* defined(__vagabond__Polymer__) */
