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

#include "Novalent.h"
#include "Bond.h"
#include "Anchor.h"
#include "Options.h"
#include "Crystal.h"
#include "Atom.h"

void Novalent::initialise()
{
	_abs = empty_vec3();
	_b = 0;
	_changedSamples = true;
}

Novalent::Novalent(AtomPtr atom) : ExplicitModel()
{
	initialise();
	_atom = atom;
	_abs = getAtom()->getInitialPosition();
	_b = getAtom()->getInitialBFactor();
}

Novalent::Novalent() : ExplicitModel()
{
	initialise();
}

Novalent::~Novalent()
{

}

std::vector<BondSample> *Novalent::getManyPositions(void *object)
{
	if (!_changedSamples)
	{
		return &_storedSamples;
	}

	CrystalPtr crystal = Options::getActiveCrystal();
	int totalPoints = crystal->getSampleNum();

	std::vector<double> occs;
	std::vector<vec3> points;
	double msq = _b / (8 * M_PI * M_PI) * 3;
	msq *= sqrt(2);
	points = ExplicitModel::makeCloud(totalPoints, msq, occs);
	_storedSamples.clear();
	
	for (int i = 0; i < points.size(); i++)
	{
		BondSample sample;
		sample.basis = make_mat3x3();
		sample.occupancy = occs[i];
		sample.torsion = 0;
		vec3 finished = _abs;
		vec3_add_to_vec3(&finished, points[i]);
		sample.start = finished;
		
		_storedSamples.push_back(sample);
	}
	
	AnchorPtr anch = getAnchor();
	if (anch && anch->motionCount() > 0)
	{
		MotionPtr mot = anch->getMotion(0);
		mot->applyTranslations(_storedSamples, true);
	}
	
	_changedSamples = false;
	return &_storedSamples;
}

AnchorPtr Novalent::getAnchor()
{
	if (_anch)
	{
		return _anch;
	}

	CrystalPtr crystal = Options::getActiveCrystal();
	std::vector<AtomPtr> chosen = crystal->getCloseAtoms(getAtom(), 5, false);

	for (size_t i = 0; i < chosen.size(); i++)
	{
		if (chosen[i]->getModel()->isBond())
		{
			_anch = ToBondPtr(chosen[i]->getModel())->getAnchor();
			return _anch;
		}
	}
	
	return AnchorPtr();
}

void Novalent::propagateChange(int depth, bool refresh)
{
	_changedSamples = true;
	ExplicitModel::propagateChange(depth, refresh);
}

std::string Novalent::shortDesc()
{
	if (_atom.expired())
	{
		return "Null_novalent";
	}
	
	return "Novalent_" + _atom.lock()->shortDesc();
}

std::string Novalent::getParserIdentifier()
{
	return shortDesc();
}

void Novalent::addProperties()
{
	addReference("atom", _atom.lock());
	addVec3Property("position", &_abs);
	addDoubleProperty("bfactor", &_b);
	
	ExplicitModel::addProperties();
}

void Novalent::postParseTidy()
{
	ExplicitModel::postParseTidy();

}

void Novalent::linkReference(BaseParserPtr object, std::string category)
{
	AtomPtr atom = ToAtomPtr(object);

	if (category == "atom")
	{
		_atom = atom;
	}
}

void Novalent::addObject(ParserPtr object, std::string category)
{

}



