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

#ifndef __cluster4x__QuickAtoms__
#define __cluster4x__QuickAtoms__

#include <map>
#include <string>
#include <vector>
#include <libsrc/vec3.h>
#include <libsrc/maths.h>
#include <libsrc/shared_ptrs.h>

typedef std::vector<vec3> Vec3Vec;
typedef std::vector<size_t> Counts;

class MtzFile;
class CAlphaView;

class QuickAtoms
{
public:
	QuickAtoms(MtzFile *file);
	
	static double compare(QuickAtoms *one, QuickAtoms *two, QuickAtoms *ave);

	void fetchAtoms();
	void addAtomsFrom(QuickAtoms *other);
	void divideThrough();
	void collapseOnTarget(vec3 target);
	
	void addSequentialAtom(std::string chain, vec3 pos);
	
	vec3 getCentre()
	{
		return _centre;
	}

	void populateCAlphaView(CAlphaView *view);
private:
	void addFromChain(QuickAtoms *other, std::string chain);
	void populatePolymer(PolymerPtr p);
	MtzFile *_file;

	std::map<std::string, Vec3Vec> _chainMap;
	std::map<std::string, Counts> _countMap;
	
	vec3 _centre;
};

#endif
