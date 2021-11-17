// vagabond
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

#ifndef __vagabond__ConfSpace__
#define __vagabond__ConfSpace__

#include "shared_ptrs.h"
#include <vector>
#include <map>

typedef std::vector<double> SpacePoint;

class ConfAxis;
class SVDBond;

class ConfSpace
{
public:
	ConfSpace(int n);
	~ConfSpace();
	
	int axisCount()
	{
		return _axes.size();
	}
	
	ConfAxis *axis(int i)
	{
		return _axes[i];
	}

	void addMolecule(MoleculePtr mol);
	void calculateFrom(PolymerPtr pol);
	void readFromFile(std::string filename);
	double totalMotionForResidue(int i);
	double weightedTorsion(int resi, bool phi, SpacePoint &sp);
private:
	void assembleContributions();
	std::map<int, std::vector<double> > _phis;
	std::map<int, std::vector<double> > _psis;

	std::vector<ConfAxis *> _axes;
	std::vector<MoleculePtr> _molecules;
	SVDBond *_svd;
};

#endif
