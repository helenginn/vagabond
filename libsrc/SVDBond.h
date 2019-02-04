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

#ifndef __vagabond__svdbond__
#define __vagabond__svdbond__

#include "FlexLocal.h"
#include "Param.h"

typedef struct
{
	Param *pWhack;
	Param *pKick;
	double *rowPtr;
} ParamSVDSet;

typedef struct
{
	double kick;
	double whack;
} BondParamPair;

typedef std::map<BondPtr, BondParamPair> BondBase;

class SVDBond
{
public:
	SVDBond(BondEffects &effects, std::vector<BondPtr> &bonds,
	        std::vector<AtomPtr> &atoms);
	
	~SVDBond();

	int numClusters()
	{
		return _num;
	}

	void addToStrategy(RefinementStrategyPtr strategy, int dir);
	void applyParameters();
	void performSVD(BondBondCC *ccs = NULL);
private:
	void prepareMatrix(double ***ptr);
	void prepareVector(double **ptr);
	void populateMatrix();
	void svdMagic();
	void report();
	void cleanupSVD(double ***ptr);
	void copyMatrix(double **from, double **to);
	
	std::vector<ParamSVDSet> _params;
	
	BondEffects _effects;
	BondBase _bondBase;
	std::vector<BondPtr> _bonds;
	std::vector<AtomPtr> _atoms;

	double **_svd;
	double **_original;
	double **_v;
	double *_w;
	
	int _num;
	double _wTotal;
};


#endif
