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

#include "Selected.h"
#include "Selected2GL.h"
#include "Monomer.h"
#include "Options.h"

Selected2GL::Selected2GL()
{
	_sType = SelectAtom;
	_conf = 0;
	_mousey = false;
	_refining = false;
}

void Selected2GL::findAtoms()
{
	ExplicitModel::useMutex();
	clearVertices();
	
	if (!_picked)
	{
		return;
	}
	
	AtomGroupPtr selected = getSelected();
	
	if (!selected || !selected->atomCount())
	{
		_shouldGetBonds = false;
		return;
	}
	
	for (int i = 0; i < selected->atomCount(); i++)
	{
		if (selected->atom(i)->getElectronCount() <= 1)
		{
			continue;
		}

		addAtom(selected->atom(i));
	}
	
	_shouldGetBonds = false;
}

bool Selected2GL::isRefinable()
{
	if (!refinableSelection()->atomCount())
	{
		return false;
	}
	
	return (_picked && _sType != SelectAtom);
}

AtomGroupPtr Selected2GL::refinableSelection()
{
	AtomGroupPtr raw = getSelected();
	
	if (!raw)
	{
		return AtomGroupPtr();
	}
	
	AtomList list = raw->findAtoms("CA");

	for (int i = 0; i < list.size(); i++)
	{
		raw->removeAtom(list[i]);
	}
	
	return raw;
}

int Selected2GL::processMolecule(MoleculePtr molecule)
{
	updateAtoms();
}

void Selected2GL::bindTextures()
{
	Vagabond2GL::bindTextures();
	bindOneTexture(pic_selected);
}

AtomGroupPtr Selected2GL::getSelected()
{
	if (!_picked)
	{
		return AtomGroupPtr();
	}

	AtomGroupPtr group = AtomGroupPtr(new AtomGroup());
	group->addAtom(_picked);

	if (_sType == SelectAtom || !_picked->getMonomer())
	{
		return group;
	}

	if (_sType == SelectResidue)
	{
		MonomerPtr monomer = _picked->getMonomer();
		if (monomer)
		{
			SidechainPtr side = monomer->getSidechain();
			return side;
		}
		
		return group;
	}
	
	if (_sType == SelectConformer)
	{
		MonomerPtr monomer = _picked->getMonomer();
		
		if (monomer)
		{
			SidechainPtr side = monomer->getSidechain();
			while (true)
			{
				AtomGroupPtr subgroup = side->subGroupForConf(_conf);

				bool heavy = false;
				for (int i = 0; i < subgroup->atomCount(); i++)
				{
					if (subgroup->atom(i)->getElectronCount() > 1)
					{
						heavy = true;
						break;
					}
				}

				if (!heavy)
				{
					_conf++;
				}
				else
				{
					return subgroup;
				}

				if (_conf >= side->conformerCount())
				{
					return group;
				}
			}
		}
		
	}
	
	return AtomGroupPtr();
}

void Selected2GL::setPicked(AtomPtr atom)
{
	if (_refining)
	{
		return;
	}

	if (!atom)
	{
		_picked = AtomPtr();
		_sType = SelectAtom;
		_shouldGetBonds = true;
		return;
	}

	AtomGroupPtr selected = getSelected();
	bool inGroup = false;

	if (_picked && _picked->getMonomer())
	{
		inGroup = _picked->getMonomer()->hasAtom(atom);
	}
	
	if (atom == _picked && _sType == SelectAtom)
	{
		_sType = SelectResidue;
		Options::statusMessage("Selected residue " 
		                       + atom->getMonomer()->getIdentifier() 
		                       + i_to_str(atom->getMonomer()->getResidueNum()));
	}
	else if (atom == _picked && _sType == SelectResidue)
	{
		if (selected->conformerCount() <= 1)
		{
			_sType = SelectAtom;
		}
		else if (selected->conformerCount() > 1)
		{
			_sType = SelectConformer;
			_conf = 0;
			
			while (_conf < selected->conformerCount() - 1 && 
			       selected->conformer(_conf) == "")
			{
				_conf++;
			}
		}
	}
	else if (inGroup && _sType == SelectConformer)
	{
		if (_conf < _picked->getMonomer()->conformerCount() - 1)
		{
			_conf++;
		}
		else
		{
			_conf = 0;
			_sType = SelectAtom;
		}
	}
	
	if (_picked != atom)
	{
		Options::statusMessage("Selected atom " + atom->shortDesc());
		_sType = SelectAtom;
		_conf = 0;
	}
	
	_picked = atom;
	_shouldGetBonds = true;
}

vec3 Selected2GL::averageModelPos()
{
	AtomGroupPtr selected = refinableSelection();

	vec3 sum = empty_vec3();
	double count = 0;
	
	for (int i = 0; i < selected->atomCount(); i++)
	{
		double last = 1;
		AtomPtr atom = selected->atom(i);
		vec3 pos = atom->getAbsolutePosition();
		vec3 model = mat4x4_mult_vec3(modelMat, pos, &last);
		
		vec3_add_to_vec3(&sum, model);
		count++;
	}

	vec3_mult(&sum, 1 / count);
	return sum;
}

void Selected2GL::focusOnGroup()
{
	AtomGroupPtr group = getSelected();
	
	if (!group)
	{
		return;
	}

	vec3 centroid = group->centroid();
	
	Options::getRuntimeOptions()->focusOnPosition(centroid);
}

void Selected2GL::manualRefine()
{
	_switch.lock();
	_refining = true;
	_switch.unlock();

	std::cout << "Starting manual refinement." << std::endl;
	CrystalPtr crystal = Options::getActiveCrystal();
	AtomGroupPtr group = refinableSelection();
	vec3 ave = averageModelPos();
	mat4x4 inv = mat4x4_inverse(modelMat);

	while (true)
	{
		_switch.lock();
		bool val = _refining;
		bool mousey = _mousey;
		vec3 ray = _ray;
		_switch.unlock();

		group->addParamType(ParamOptionMaxTries, 1.0);
		group->addParamType(ParamOptionNumBonds, 4.0);
		group->addParamType(ParamOptionTorsion, 0.5);
		
		if (!val)
		{
			std::cout << "Ending manual refinement." << std::endl;
			break;
		}
		
		if (mousey)
		{
			group->addParamType(ParamOptionTorsion, 0.2);
			double z_ratio = ave.z / ray.z;
			vec3_mult(&ray, z_ratio);
			ray = mat4x4_mult_vec(inv, ray);
			double max_weight = 0;
			
			for (int i = 0; i < group->atomCount(); i++)
			{
				AtomPtr atom = group->atom(i);
				vec3 pos = atom->getAbsolutePosition();
				vec3 model = mat4x4_mult_vec(modelMat, pos);
				vec3 diff = vec3_subtract_vec3(model, ray);
				double weight = 1 / vec3_length(diff);
				atom->setTargetPosition(ray, weight);
				
				if (max_weight < weight)
				{
					max_weight = weight;
				}
			}

			for (int i = 0; i < group->atomCount(); i++)
			{
				AtomPtr atom = group->atom(i);
				double w = atom->getTargetWeight();
			
				if (w < max_weight * 0.2)	
				{
					atom->setTargetPosition(ray, 0);
				}
			}

			group->refine(crystal, RefinementMouse);
		}
		else
		{
			group->refine(crystal, RefinementFine);
		}
	}

	group->clearParams();

}

void Selected2GL::setMouseRay(vec3 ray)
{
	_switch.lock();
	_ray = ray;
	_switch.unlock();
}

void Selected2GL::setMouseRefinement(bool val)
{
	_switch.lock();
	_mousey = val;
	_switch.unlock();
}

void Selected2GL::splitSelected()
{
	if (!_picked)
	{
		Options::statusMessage("Split failed (no selection)");
		return;
	}
	if (_sType == SelectAtom)
	{
		Options::statusMessage("Split failed (only one atom selected).");
		return;
	}

	AtomGroupPtr group = getSelected();
	
	if (group->conformerCount() == 1)
	{
		/* Acceptable */
		
		AtomPtr cb = group->findAtom("CB");
		
		if (!cb)
		{
			Options::statusMessage("Split failed (no Cb atom).");
			return;
		}
		
		ModelPtr model = cb->getModel();
		
		if (!model->isBond())
		{
			Options::statusMessage("Split failed (Cb is not bonded).");
			return;
		}
		
		ToBondPtr(model)->splitBond();
		Options::statusMessage("Split success (" + cb->shortDesc() + ").");
	}
	else
	{
		Options::statusMessage("Split failed (too many conformers selected.)");
	}
}
