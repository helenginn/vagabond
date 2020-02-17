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
#include "../libsrc/Monomer.h"
#include "../libsrc/Options.h"
#include "../libsrc/Polymer.h"

Selected2GL::Selected2GL()
{
	_sType = SelectAtom;
	_conf = 0;
	_mousey = false;
	_refining = false;
	_kicking = false;
	_adding = false;
}

void Selected2GL::findAtoms()
{
	ExplicitModel::useMutex();
	clearVertices();
	
	if (!getPicked())
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

AtomGroupPtr Selected2GL::refinableSelection()
{
	AtomGroupPtr raw = getSelected();
	
	if (!raw || _sType == SelectAtom)
	{
		return AtomGroupPtr();
	}
	
	AtomList list = raw->findAtoms("CA");
	AtomList more = raw->findAtoms("HA");
	list.reserve(list.size() + more.size());
	list.insert(list.end(), more.begin(), more.end());

	if (_sType == SelectMonomer || _sType == SelectMonConf)
	{
		return raw;
	}

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
	if (!getPicked())
	{
		return AtomGroupPtr();
	}

	AtomGroupPtr group = AtomGroupPtr(new AtomGroup());
	
	for (int i = 0; i < _picked.size(); i++)
	{
		AtomPtr pick = _picked[i];
		
		if (!pick)
		{
			continue;
		}
		
		if (_sType == SelectAtom)
		{
			group->addAtom(pick);

			if (pick->isHeteroAtom())
			{
				continue;
			}
		}

		if (_sType == SelectMonomer)
		{
			MonomerPtr monomer = pick->getMonomer();
			group->addAtomsFrom(monomer);
		}

		if (_sType == SelectMonConf)
		{
			MonomerPtr monomer = pick->getMonomer();
			AtomGroupPtr subgroup = monomer->subGroupForConf(_conf);
			group->addAtomsFrom(subgroup);
		}

		if (_sType == SelectSidechain)
		{
			MonomerPtr monomer = pick->getMonomer();
			if (monomer)
			{
				SidechainPtr side = monomer->getSidechain();
				group->addAtomsFrom(side);
			}
		}

		if (_sType == SelectConformer)
		{
			MonomerPtr monomer = pick->getMonomer();

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
						group->addAtomsFrom(subgroup);
						break;
					}

					if (_conf >= side->conformerCount())
					{
						break;
					}
				}
			}

		}
	}
	
	return group;
}

void Selected2GL::addPicked(AtomPtr atom, bool preserveType)
{
	/* If we didn't select an atom, but we still had shift pressed,
	 * probably a mistake, so ignore entirely */
	if (!atom)
	{
		return;
	}

	bool inGroup = false;
	_shouldGetBonds = true;
	AtomGroupPtr selected = getSelected();

	/* If we already have a selection and want to know if it's in the
	 * current group or not */
	if (selected)
	{
		inGroup = selected->hasAtom(atom);
	}
	
	if (inGroup && _sType == SelectMonomer)
	{
		std::cout << "In group." << std::endl;
		_conf = 0;
		
		if (_conf == 0)
		{
			_conf++;
		}
		
		if (_conf >= selected->conformerCount())
		{
			_sType = SelectMonomer;
			_conf = 0;
			return;
		}
		
		if (selected->conformer(_conf) == "a" ||
		    selected->conformer(_conf) == "A")
		{
			_conf++;
		}
		
		_sType = SelectMonConf;
		Options::statusMessage("Selected conformer "
		                       + selected->conformer(_conf), false);

		return;
	}
	else if (inGroup && _sType == SelectMonConf)
	{
		_conf++;
		
		if (_conf >= selected->conformerCount())
		{
			_sType = SelectMonomer;
			_conf = 0;
			return;
		}
		
		
		Options::statusMessage("Selected conformer "
		                       + selected->conformer(_conf), false);
	}
	else
	{
		_picked.push_back(atom);

		if (!preserveType && !atom->isHeteroAtom())
		{
			_sType = SelectMonomer;
		}
	}
}

void Selected2GL::setPicked(AtomPtr atom, bool preserveType)
{
	/* Don't change selections when refining */
	if (_refining)
	{
		return;
	}

	/* If we have pressed shift we deal with multi-atom selections */
	if (_adding)
	{
		addPicked(atom, preserveType);
		
		return;
	}
	
	_shouldGetBonds = true;

	/* If we have picked empty, we are not pressing shift anymore
	 * but still have a multi-atom selection */
	if (hasMultiAtomSelection() && !atom)
	{
		_picked.clear();
		_picked.push_back(AtomPtr());
		return;
	}
	/* If we have picked a real atom, we are not pressing shift anymore
	 * but still have a multi-atom selection */
	else if (hasMultiAtomSelection() && atom)
	{
		AtomPtr first = getPicked();
		_picked.clear();
		_picked.push_back(first);
	}
	
	/* Now we don't have a multi-atom selection */

	/* We picked empty */
	if (!atom)
	{
		_sType = SelectAtom;

		if (_picked.size())
		{
			_picked[0] = atom;
		}
		else
		{
			_picked.push_back(atom);
		}
		return;
	}

	AtomGroupPtr selected = getSelected();
	bool inGroup = false;

	if (getPicked() && !getPicked()->isHeteroAtom())
	{
		MonomerPtr mon = getPicked()->getMonomer();
		
		if (mon)
		{
			inGroup = getPicked()->getMonomer()->hasAtom(atom);
		}
	}
	
	if (atom == getPicked() && _sType == SelectAtom
	    && !getPicked()->isHeteroAtom())
	{
		_sType = SelectSidechain;
		Options::statusMessage("Selected residue " 
		                       + atom->getMonomer()->getIdentifier() 
		                       + i_to_str(atom->getMonomer()->getResidueNum())
		                       + " from chain " + 
		                       atom->getMolecule()->getChainID(), false);
	}
	else if (inGroup && _sType == SelectSidechain && !preserveType)
	{
		if (selected->conformerCount() <= 1)
		{
			_sType = SelectAtom;
		}
		else if (selected->conformerCount() > 1)
		{
			_sType = SelectConformer;
			_conf = 0;
			
			if (atom->getMonomer()->conformer(_conf) == "")
			{
				_conf++;
			}
			
			std::string status = "Selected residue " 
			+ atom->getMonomer()->getIdentifier() 
			+ i_to_str(atom->getMonomer()->getResidueNum()) + "_"
			+ atom->getMonomer()->conformer(_conf);
			
			Options::statusMessage(status, false);
		}
	}
	else if (inGroup && _sType == SelectConformer && !preserveType)
	{
		_conf++;

		if (_conf >= atom->getMonomer()->conformerCount())
		{
			_sType = SelectSidechain;
			_conf = 0;
		}
		else
		{
			std::string status = "Selected residue " 
			+ atom->getMonomer()->getIdentifier() 
			+ i_to_str(atom->getMonomer()->getResidueNum()) + "_"
			+ atom->getMonomer()->conformer(_conf);
		}
	}
	
	if (getPicked() != atom)
	{
		Options::statusMessage("Selected atom " + atom->shortDesc(), false);
		
		if (!preserveType)
		{
			_sType = SelectAtom;
		}

		_conf = 0;
	}
	
	/* Change the only atom if one exists already */
	if (_picked.size())
	{
		_picked[0] = atom;
	}
	else
	{
		_picked.push_back(atom);
	}
}

vec3 Selected2GL::averageModelPos()
{
	AtomGroupPtr selected = refinableSelection();
	
	if (!selected)
	{
		return empty_vec3();
	}

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

void Selected2GL::toggleKicks()
{
	_switch.lock();
	_kicking = !_kicking;
	
	if (_kicking)
	{
		Options::statusMessage("Add kicking to refinement.", false);
	}
	else
	{
		Options::statusMessage("Remove kicking in refinement.", false);
	}

	_switch.unlock();
}

void Selected2GL::advanceMonomer(int dir)
{
	int res = 0;
	PolymerPtr pol;

	if (getPicked() && 
	    getPicked()->getMolecule() && 
	    getPicked()->getMolecule()->isPolymer())
	{
		pol = ToPolymerPtr(getPicked()->getMolecule());
		res = getPicked()->getResidueNum();
		res += dir;
	
	}
	else
	{
		CrystalPtr crystal = Options::getActiveCrystal();
		
		for (int i = 0; i < crystal->moleculeCount(); i++)
		{
			if (crystal->molecule(i)->isPolymer())
			{
				pol = ToPolymerPtr(crystal->molecule(i));
				res = pol->monomerBegin();
			}
		}
	}

	MonomerPtr mon = pol->getMonomer(res);

	if (mon)
	{
		AtomPtr atom = mon->findAtom("CA");
		setPicked(atom, true);
	}

	focusOnGroup();
}

void Selected2GL::cancelRefine()
{
	Options::statusMessage("Cancelling refinement...", false);
	_switch.lock();
//	_kicking = false;
	_refining = false;
	_switch.unlock();
}

void Selected2GL::manualRefine()
{
	if (!getPicked())
	{
		std::cout << "No atoms selected for manual refinement." << std::endl;
		return;
	}

	MoleculePtr mol = getPicked()->getMolecule();
	if (!mol->isPolymer())
	{
		std::cout << "Need to refine atoms from polymer" << std::endl;
		return;
	}

	_switch.lock();
	_refining = true;
	_switch.unlock();

	Options::statusMessage("Starting manual refinement.", false);
	CrystalPtr crystal = Options::getActiveCrystal();
	AtomGroupPtr group = refinableSelection();
	
	if (!group)
	{
		return;
	}
	
	bool terminal = (group->beyondGroupAtoms().size() == 0);
	int begin = 0; int end = 0;
	
	if (!terminal)
	{
		group->boundingMonomers(&begin, &end);
	}

	/* These will take on old values for comparison */
	bool refining = false;
	bool mousey = false;
	bool kicking = false;

	while (true)
	{
		bool was_mousey = mousey;

		_switch.lock();
		refining = _refining;
		mousey = _mousey;
		kicking = _kicking;
		vec3 ray = _ray;
		_switch.unlock();

		vec3 ave = averageModelPos();
		mat4x4 inv = mat4x4_inverse(modelMat);

		group->clearParams();
		group->addParamType(ParamOptionMaxTries, 1.0);
		group->addParamType(ParamOptionNumBonds, 4.0);
		group->addParamType(ParamOptionTorsion, 0.5);
		
		if (kicking && !mousey)
		{
			group->addParamType(ParamOptionKick, 0.2);
		}
		
		if (!refining)
		{
			Options::statusMessage("Ending manual refinement.", false);
			break;
		}
		
		if (mousey && !was_mousey)
		{
			Options::statusMessage("Mouse-driven refinement.", false);
		}
		else if (!mousey && was_mousey)
		{
			Options::statusMessage("Refining to electron density.", false);
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
		else if (terminal)
		{
			group->addParamType(ParamOptionBondAngle, 0.5);
			group->refine(crystal, RefinementFine);
			bool changed = group->didChange();

			group->addParamType(ParamOptionMaxTries, 1.0);
			group->addParamType(ParamOptionNumBonds, 0.0);
			group->addParamType(ParamOptionOccupancy, 1.0);

			group->refine(crystal, RefinementFine);
			changed &= group->didChange();

			if (!changed && false)
			{
				std::cout << "No change from previous refinement. Done!"
				<< std::endl;
				Options::statusMessage("No more improvement, stopping.", false);
				cancelRefine();	
				break;
			}
		}
		else if (!terminal)
		{
			PolymerPtr pol = ToPolymerPtr(mol);
			pol->refineFromFarRegion(begin, end, crystal);
			
			if (!pol->didChange())
			{
				std::cout << "No change from previous refinement."
				<< std::endl;
			}
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

	if (val)
	{
		Options::statusMessage("Preparing mouse refinement...", false);
	}
	else
	{
		Options::statusMessage("Refining to electron density.", false);
	}
}

void Selected2GL::splitSelected()
{
	if (!getPicked())
	{
		Options::statusMessage("Split failed (no selection)", false);
		return;
	}
	if (_sType == SelectAtom)
	{
		Options::statusMessage("Split failed (only one atom selected).", false);
		return;
	}

	AtomGroupPtr group = getSelected();
	
	if (_sType == SelectSidechain || _sType == SelectConformer)
	{
		if (group->conformerCount() == 1)
		{
			/* Acceptable */

			AtomPtr cb = group->findAtom("CB");

			if (!cb)
			{
				Options::statusMessage("Split failed (no Cb atom).", false);
				return;
			}

			ModelPtr model = cb->getModel();

			if (!model->isBond())
			{
				Options::statusMessage("Split failed (Cb is not bonded).", false);
				return;
			}

			ToBondPtr(model)->splitBond();
			Options::statusMessage("Split success (" + cb->shortDesc() + ").", false);
		}
		else
		{
			Options::statusMessage("Split failed (too many conformers "
			                       "selected.)", false);
		}
		
		return;
	}
	
	if (_sType != SelectMonomer)
	{
		Options::statusMessage("Split failed (unknown selection option.)", false);
		return;
	}
	
	AtomList blockable = group->beyondGroupAtoms();
	
	for (int i = 0; i < blockable.size(); i++)
	{
		std::cout << "Split blocking " << blockable[i] << std::endl;
		ToBondPtr(blockable[i]->getModel())->setSplitBlock();
	}
	
	AtomList tops = group->topLevelAtoms();
	
	for (int i = 0; i < tops.size(); i++)
	{
		if (tops[i]->getModel()->isBond())
		{
			ToBondPtr(tops[i]->getModel())->splitBond();
		}
	}

	Options::statusMessage("Split success (" + tops[0]->shortDesc() + ").", false);
}

void Selected2GL::deleteSelected()
{
	if (!getPicked())
	{
		return;
	}
	
	if (getPicked()->isHeteroAtom())
	{
		Options::statusMessage("Deleted " + getPicked()->shortDesc(), false);
		Options::getActiveCrystal()->removeAtom(getPicked());

		_picked.clear();
	}
}


void Selected2GL::selectResidue(std::string chain, int resNum)
{
	CrystalPtr crystal = Options::getActiveCrystal();
	
	if (!chain.length())
	{
		if (getPicked())
		{
			chain = getPicked()->getMolecule()->getChainID();
		}
		else
		{
			if (crystal->atomCount())
			{
				chain = crystal->atom(0)->getMolecule()->getChainID();
			}
			else
			{
				return;
			}
		}
	}

	for (int i = 0; i < crystal->moleculeCount(); i++)
	{
		MoleculePtr molecule = crystal->molecule(i);

		if (!molecule->isPolymer())
		{
			continue;
		}
		
		if (chain.length() > 0 && molecule->getChainID()[0] != chain[0])
		{
			continue;
		}
		
		PolymerPtr pol = ToPolymerPtr(molecule);
		
		if (resNum == -INT_MAX)
		{
			resNum = pol->monomerBegin();
		}

		if (pol->hasResidue(resNum))
		{
			MonomerPtr mon = pol->getMonomer(resNum);

			if (mon)
			{
				AtomPtr atom = mon->findAtom("CA");
				setPicked(atom, true);
			}

			focusOnGroup();
			return;
		}
	}
	
	std::cout << "Did not find residue " << chain << " " << resNum << std::endl;
}
