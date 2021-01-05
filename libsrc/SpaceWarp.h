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

#include "shared_ptrs.h"
#include <hcsrc/vec3.h>
#include <hcsrc/Any.h>
#include <map>

typedef std::vector<BondPtr> BondList;

class FlexGlobal;

class SpaceWarp
{
public:
	SpaceWarp(VagFFTPtr fft);
	~SpaceWarp();
	
	void recalculate(VagFFTPtr data);
	
	double evaluate();
	
	static double getScore(void *object)
	{
		return static_cast<SpaceWarp *>(object)->evaluate();
	}
	
	void addBond(BondPtr bond)
	{
		_bonds.push_back(bond);
	}
	
	void addRefinedAtom(AtomPtr atom);

	VagFFTPtr getFFT()
	{
		return _fft;
	}
	
	void setCalculated(VagFFTPtr);

	vec3 getWarp(int n)
	{
		return _warps[n];
	}
	
	void setActiveAtoms(AtomGroupPtr atoms);

	static double staticScore(void *object)
	{
		return static_cast<SpaceWarp *>(object)->score();
	}
	double score();
	void svd();
private:
	vec3 bondEffect(vec3 pos, BondPtr b);
	void setupSVD();
	void wipeSVD();
	void populateSVD(AtomPtr atom, BondPtr bond, int interaction);
	void populateSVD(AtomPtr a);
	void cleanupSVD();
	void addTargets();
	void addBondAtoms(BondPtr bond);
	void initAtom(AtomPtr atom);
	void nudgeAtom(AtomPtr atom);
	void nudge();
	vec3 target(AtomPtr atom);
	double cost();
	VagFFTPtr _fft;
	VagFFTPtr _data;
	
	std::map<AtomPtr, std::map<BondPtr, int> > _interactions;
	
	AtomGroupPtr _atomsForSVD;
	AtomGroupPtr _atoms;
	AtomGroupPtr _all;
	std::vector<BondPtr> _bonds;
	
	/* offsets in voxels */
	vec3 *_warps;
	
	std::vector<AnyPtr> _varying;
	
	vec3 _activePos;
	
	FlexGlobal *_target;

	double _magnitude;
	double *_matrix;
	double *_w;
	double *_v;
	double *_t;
	double *_weights;
	double *_torsions;

	double **_matPtrs;
	double **_vPtrs;
	
	int _atomCounter;
	int _bondCounter;
};

