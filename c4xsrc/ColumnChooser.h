// cluster4x
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

#ifndef __cluster4x__ColumnChooser__
#define __cluster4x__ColumnChooser__

#include <vector>
#include <map>
#include "MtzFFTPtr.h"

class Group;

class ColumnChooser
{
public:
	ColumnChooser();
	
	void addTargetGroup(Group *g);
	
	double evaluate(int enabled, bool work = false);
	void prune();
private:
	void enableAllColumns();
	bool pruneCycle();
	void randomFree(double prop = 0.4);
	void loadSavedSet();
	void saveColumnSet();

	std::vector<Group *> _groups;
	std::vector<MtzFFTPtr> _mtzs;
	std::map<MtzFFTPtr, bool> _freeMap;
	std::map<int, bool> _saved;

};

#endif

