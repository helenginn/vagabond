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
	Param *pPhi;
	Param *pTorsion;
	double *rowPtr;
} ParamSVDSet;

typedef struct
{
	double kick; /* Baseline kick, minus adjustments */
	double whack; /* Baseline whack, minus adjustments */
	double torsion; /* Baseline torsion, minus adjustments */
	double phi; /* Baseline phi, minus adjustments */
} BondParamPair;

typedef std::map<BondPtr, BondParamPair> BondBase;

class SVDBond
{
public:
	SVDBond(std::vector<BondPtr> &bonds, std::vector<AtomPtr> &atoms);
	
	~SVDBond();

	int numClusters()
	{
		return _num;
	}
	
	void setSilent(bool s)
	{
		_silent = s;
	}
	
	void setDoTorsion(bool t)
	{
		_doTorsion = t;
	}
	
	void setDoAngles(bool a)
	{
		_doAngles = a;
	}

	void addToStrategy(RefinementStrategyPtr strategy, double mult,
	                   bool phi = false);
	void applyParameters();
	void performSVD();
private:
	void prepareMatrix(double ***ptr);
	void prepareVector(double **ptr);
	void svdMagic();
	void report();
	void cleanupSVD(double ***ptr);
	void copyMatrix(double **from, double **to);
	void compareBonds();
	double compareTorsionWithAngle(BondPtr a, BondPtr b);
	double compareTorsions(BondPtr a, BondPtr b);
	double compareForKicks(BondPtr a, BondPtr b);
	double compareKicks(BondPtr a, BondPtr b);
	void writeMatrix();
	void determineInteractions();
	void interactionsForBond(BondPtr bond);
	
	std::vector<ParamSVDSet> _params;
	
	/* basic parameters before refinement additions */
	BondBase _bondBase;

	std::map<AtomPtr, std::map<BondPtr, int> > _interactions;

	std::vector<BondPtr> _bonds;
	std::vector<AtomPtr> _atoms;

	double **_svd;
	double **_original;
	double **_v;
	double *_w;
	
	int _num;
	double _wTotal;
	bool _doTorsion;
	bool _doAngles;
	bool _silent;
};

vec3 bond_effect_on_pos(vec3 atom_pos, mat3x3 &bond_basis, vec3 &bond_pos);

#endif
