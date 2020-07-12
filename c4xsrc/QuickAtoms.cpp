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

#include "QuickAtoms.h"
#include "MtzFile.h"
#include "CAlphaView.h"
#include <libsrc/Crystal.h>
#include <libsrc/Atom.h>
#include <libsrc/Monomer.h>
#include <libsrc/Polymer.h>
#include <libsrc/SymMate.h>

QuickAtoms::QuickAtoms(MtzFile *file)
{
	_file = file;
}

void QuickAtoms::addFromChain(QuickAtoms *other, std::string chain)
{
	Vec3Vec theirs = other->_chainMap[chain];
	Vec3Vec mine = _chainMap[chain];
	Counts counts = _countMap[chain];
	
	if (theirs.size() > mine.size())
	{
		mine.resize(theirs.size());
		counts.resize(theirs.size());
	}
	
	for (size_t i = 0; i < theirs.size(); i++)
	{
		if (theirs[i].x != theirs[i].x)
		{
			continue;
		}

		vec3_add_to_vec3(&mine[i], theirs[i]);
		counts[i] += 1;
	}
	
	_chainMap[chain] = mine;
	_countMap[chain] = counts;
}

void QuickAtoms::addAtomsFrom(QuickAtoms *other)
{
	if (other->_chainMap.size() == 0)
	{
		return;
	}

	for (std::map<std::string, Vec3Vec>::iterator it
	     = other->_chainMap.begin();
	     it != other->_chainMap.end(); it++)
	{
		std::string chain = it->first;
		addFromChain(other, chain);
	}
}

void QuickAtoms::divideThrough()
{
	_centre = empty_vec3();
	int total = 0;

	for (std::map<std::string, Vec3Vec>::iterator it
	     = _chainMap.begin();
	     it != _chainMap.end(); it++)
	{
		std::string chain = it->first;
		
		for (size_t i = 0; i < _chainMap[chain].size(); i++)
		{
			double scale = 1 / (double)_countMap[chain][i];
			
			vec3_mult(&_chainMap[chain][i], scale);
			_countMap[chain][i] = 1;
			vec3 v = _chainMap[chain][i];
			
			if (v.x == v.x)
			{
				vec3_add_to_vec3(&_centre, _chainMap[chain][i]);
				total++;
			}
		}
	}
	
	vec3_mult(&_centre, 1 / (double)total);
}

void QuickAtoms::fetchAtoms()
{
	if (_file == NULL)
	{
		return;
	}

	if (_chainMap.size() > 0)
	{
		_chainMap.clear();
		_countMap.clear();
	}

	CrystalPtr c = _file->getCrystal();
	
	if (!c)
	{
		return;
	}

	for (size_t j = 0; j < c->moleculeCount(); j++)
	{
		MoleculePtr mol = c->molecule(j);

		if (!mol->isPolymer())
		{
			continue;
		}

		PolymerPtr p = ToPolymerPtr(mol);

		populatePolymer(p);
	}
	
	_centre = c->centroid();
}

void QuickAtoms::collapseOnTarget(vec3 target)
{
	CrystalPtr c = _file->getCrystal();
	
	if (!c)
	{
		return;
	}
	
	SymMate mate(c);
	mate.findSymop(target);
	mate.applySymops(c);

	fetchAtoms();
}

void QuickAtoms::populatePolymer(PolymerPtr p)
{
	Vec3Vec positions;
	std::string myChain;
	myChain = p->getChainID()[0];
	
	if (_chainMap.count(myChain))
	{
		positions = _chainMap[myChain];
	}

	size_t ending = positions.size();
	positions.resize(p->monomerEnd());
	
	for (size_t i = ending; i < positions.size(); i++)
	{
		/* default nan */
		vec3_mult(&positions[i], std::nan(""));
	}

	for (int i = p->monomerBegin(); i < p->monomerEnd(); i++)
	{
		if (i < 0)
		{
			std::cout << "Warning: ignoring negative residue "
			<< i << " from chain " << p->getChainID() << std::endl;
			continue;
		}

		MonomerPtr m = p->getMonomer(i);
		
		if (!m)
		{
			continue;
		}
		
		AtomPtr a = m->findAtom("CA");

		if (!a)
		{
			continue;
		}

		int id = a->getResidueNum();

		vec3 pos = a->getPDBPosition();
		positions[id] = pos;
	}
	
	_chainMap[myChain] = positions;
}

double QuickAtoms::compare(QuickAtoms *one, QuickAtoms *two, QuickAtoms *ave)
{
	CorrelData cd = empty_CD();

	for (std::map<std::string, Vec3Vec>::iterator it
	     = ave->_chainMap.begin();
	     it != ave->_chainMap.end(); it++)
	{
		std::string chain = it->first;
		
		if (one->_chainMap.count(chain) == 0 ||
		    two->_chainMap.count(chain) == 0)
		{
			continue;
		}
		
		for (size_t i = 0; i < ave->_chainMap[chain].size(); i++)
		{
			if (i >= one->_chainMap[chain].size() ||
			    i >= two->_chainMap[chain].size())
			{
				continue;
			}

			vec3 avePos = ave->_chainMap[chain][i];
			vec3 onePos = one->_chainMap[chain][i];
			vec3 twoPos = two->_chainMap[chain][i];
			
			vec3_subtract_from_vec3(&onePos, avePos);
			vec3_subtract_from_vec3(&twoPos, avePos);

			add_to_CD(&cd, onePos.x, twoPos.x);
			add_to_CD(&cd, onePos.y, twoPos.y);
			add_to_CD(&cd, onePos.z, twoPos.z);
		}
	}

	double cc = evaluate_CD(cd);
	if (cc < 0)
	{
		cc = 0;
	}

	return cc;
}

void QuickAtoms::populateCAlphaView(CAlphaView *view)
{
	vec3 nanVec = make_vec3(std::nan(""), std::nan(""), std::nan(""));

	for (std::map<std::string, Vec3Vec>::iterator it
	     = _chainMap.begin();
	     it != _chainMap.end(); it++)
	{
		std::string chain = it->first;
		
		for (size_t i = 0; i < _chainMap[chain].size(); i++)
		{
			vec3 pos = _chainMap[chain][i];
			view->addCAlpha(pos);
		}

		view->addCAlpha(nanVec);
	}

	view->addCAlpha(nanVec);
}
