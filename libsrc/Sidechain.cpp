//
//  Sidechain.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Shouter.h"
#include "Sidechain.h"
#include "Sampler.h"
#include "Bond.h"
#include "Atom.h"
#include <iostream>
#include "Monomer.h"
#include "FileReader.h"
#include "Absolute.h"
#include "Crystal.h"
#include "Backbone.h"
#include "../libinfo/RotamerTable.h"

void Sidechain::refine(CrystalPtr target, RefinementType rType)
{
	if (!canRefine()) return;
	
	bool lowRes = (target->getMaxResolution() > 2.5);
	
	if (!paramCount())
	{
		double range = 2.;

		if (_timesRefined >= 5)
		{
			range = 0.2;
		}

		switch (rType)
		{
			case RefinementModelPos:
			addParamType(ParamOptionTorsion, range);
			addParamType(ParamOptionNumBonds, 3);
			addParamType(ParamOptionMaxTries, 60);
			break;

			case RefinementSavedPos:
			addParamType(ParamOptionTorsion, range);
			addParamType(ParamOptionNumBonds, 3);
			break;

			case RefinementSidePos:
			addParamType(ParamOptionTorsion, 0.1);
			addParamType(ParamOptionNumBonds, 3);
			break;

			case RefinementSidechain:
			addParamType(ParamOptionTorsion, 0.1);
			addParamType(ParamOptionKick, 0.5);
			
			if (!lowRes)
			{
				addParamType(ParamOptionMagicAngles, 30.0);
			}

			addParamType(ParamOptionNumBonds, 3);
			break;

			default:
			break;
		}
		
		if (_timesRefined >= 5)
		{
			switch (rType)
			{
				case RefinementModelPos:
				addParamType(ParamOptionTorsion, 5);
				addParamType(ParamOptionSVD, 1);
				addParamType(ParamOptionCycles, 30);
				addParamType(ParamOptionMaxTries, 10);
				break;

				case RefinementSavedPos:
				addParamType(ParamOptionBondAngle, range / 1);
				break;

				case RefinementSidechain:
				if (!lowRes)
				{
//					addParamType(ParamOptionBondAngle, 0.1);
				}

				default:
				break;
			}
		}
	}

	if (rType == RefinementSidechain || rType == RefinementSidePos)
	{
		std::cout << getMonomer()->getResCode() << std::flush;
		rType = RefinementFine;
	}

	
	MonomerPtr monomer = getMonomer();
	BackbonePtr backbone = monomer->getBackbone();

	if (backbone && rType == RefinementFine)
	{
		addIncludeForRefinement(backbone);
	}

	AtomGroup::refine(target, rType);
	clearParams();
	clearIncludeForRefinements();
}

void Sidechain::splitConformers(int count)
{
	if (count < 0)
	{
		count = conformerCount();
	}
	
	bool new_conf = (count > conformerCount());

	if (count <= 1) return;

	AtomPtr start = findAtom("CB");
	std::string duplStart = "CB";

	if (!start || !start->getModel()->isBond())
	{
		return;
	}

	AtomList caAtoms = findAtoms("CA");

	BondPtr blockBond;
	
	if (caAtoms.size() >= 1)
	{
		MonomerPtr monomer = getMonomer();
		AtomPtr blocker;
		
		if (monomer->isAfterAnchor())
		{
			blocker = monomer->findAtom("C");
		}
		else
		{
			blocker = monomer->findAtom("N");
		}
		
		if (blocker)
		{
			ModelPtr model = blocker->getModel();
			
			if (model->isBond())
			{
				blockBond = ToBondPtr(model);
			}
			else
			{
				std::cout << "Warning! Atom " << blocker->shortDesc();
				std::cout << " is not bonded." << std::endl;
			}
		}
		else
		{
				std::cout << "Warning! No blocker for ";
				std::cout << getMonomer()->getIdentifier() << std::endl;
		}

		start = findAtom("CA");
		duplStart = "CA";
	}

	BondPtr bond = ToBondPtr(start->getModel());

	for (int i = 1; i < count; i++)
	{
		if (blockBond)
		{
			blockBond->setSplitBlock();
		}

		bond->splitBond();
	}

	AtomList atoms = findAtoms(duplStart);

	for (int i = 0; i < atoms.size(); i++)
	{
		if (new_conf)
		{
			break;
		}
		
		AtomPtr atom = atoms[i];
		if (atom->getModel()->isBond())
		{
			BondPtr bond = ToBondPtr(atom->getModel());
			double origOcc = atom->getOriginalOccupancy();
			Bond::setOccupancy(&*bond, origOcc);
		}
	}

	for (int i = 0; i < getMonomer()->atomCount(); i++)
	{
		AtomPtr atom = getMonomer()->atom(i);

		if (atom->getModel()->isAbsolute() && atom->getAlternativeConformer().length())
		{
			atom->setWeighting(0);
		}
	}
}

void Sidechain::parameteriseAsRotamers()
{
	RotamerTable *table = RotamerTable::getTable();
	std::string res = getMonomer()->getIdentifier();

	std::vector<Rotamer> rotamers = table->rotamersFor(res);

	if (!rotamers.size())
	{
		return;
	}

	// do the stuff
	_rotamerised = true;
	_canRefine = false;

	splitConformers(rotamers.size());

	// rotamers.size() should equal size of atoms list.
	std::cout << "Multi-rotamer model for residue " << i_to_str(_resNum) << std::endl;

	for (int i = 0; i < rotamers.size(); i++)
	{
		for (int j = 0; j < rotamers[i].torsions.size(); j++)
		{
			TorsionAngle angle = rotamers[i].torsions[j];
			double torsionValue = deg2rad(angle.torsion);
			std::string atomID = angle.secondAtom;

			AtomList atoms = findAtoms(atomID);

			if (atoms.size() != rotamers.size())
			{
				std::cout << "Ooooo no" << std::endl;
			}

			AtomPtr atom = atoms[i];
			if (!atom->getModel()->isBond())
			{
				continue;
			}

			BondPtr bond = ToBondPtr(atom->getModel());

			Bond::setTorsion(&*bond, torsionValue);

			if (bond->getMinor()->getAtomName() == "CB")
			{
				BackbonePtr bb = getMonomer()->getBackbone();

				if (!bb)
				{
					std::cout << "Cannot find backbone for " << getMonomer()->getResidueNum() << std::endl;
					continue;
				}

				AtomPtr heavy = bb->findAtom("N");

				if (!heavy)
				{
					std::cout << "Cannot find nitrogen atom of backbone." << std::endl;
					continue;
				}

				bond->recalculateTorsion(heavy, torsionValue);
			}

			bond->setFixed(true);
		}
	}

	AtomList cbs = findAtoms("CB");

	for (int i = 0; i < cbs.size(); i++)
	{
		AtomPtr cb = cbs[i];
		if (!cb->getModel()->isBond()) continue;

		BondPtr bond = ToBondPtr(cb->getModel());
		bond->setFixed(false);
		bond->setUsingTorsion(false);

		double occ = rotamers[i].allOccupancy;
		Bond::setOccupancy(&*bond, occ);
	}
}

void Sidechain::addProperties()
{
	addBoolProperty("rotamerised", &_rotamerised);
	addBoolProperty("can_refine", &_canRefine);
	addIntProperty("res_num", &_resNum);
	AtomGroup::addProperties();
}
