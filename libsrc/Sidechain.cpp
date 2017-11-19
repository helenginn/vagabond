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
#include "Backbone.h"
#include "../libinfo/RotamerTable.h"

bool Sidechain::shouldRefineMagicAxis(BondPtr bond)
{
	return (bond->getMinor()->getAtomName() == "CB");
}

void Sidechain::refine(CrystalPtr target, RefinementType rType)
{
	if (!_rotamerised)
	{
		if (rType != RefinementFine)
		{
			AtomGroup::refine(target, rType);
		}

		return;
	}



	if (!canRefine()) return;

	if (rType != RefinementFine)
	{
		return;
	}


	AtomList atoms = findAtoms("CB");

	for (int i = 0; i < atoms.size(); i++)
	{
		setupNelderMead();
		setCrystal(target);
		setScoreType(ScoreTypeMultiply);

		if (atoms[i].expired())
		{
			continue;
		}

		AtomPtr atom = atoms[i].lock();
		std::string conformer = atom->getAlternativeConformer();

		addSampledAtoms(shared_from_this(), conformer);
		double conformerScore = score(this);

		std::cout << "conformer " << conformer <<
		" " << conformerScore << std::endl;
	}

	/*
	setScoreType(ScoreTypeRFactor);
	addRotamer(this, 4.0, 0.05);

	//	 setSilent();

	setJobName("rotamer_" + getMonomer()->getResCode()
			   + i_to_str(getMonomer()->getResidueNum()));

	sample();
	 */
}

void Sidechain::fixBackboneTorsions(AtomPtr betaTorsion)
{
	AtomPtr atom = findAtom("CB");

	if (!atom)
	{
		return;
	}

	ModelPtr model = atom->getModel();

	if (model->isBond())
	{
		ToBondPtr(model)->setTorsionAtoms(betaTorsion);
	}
}

void Sidechain::setInitialDampening()
{
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->isBackbone())
		{
			continue;
		}

		BondPtr bond = BondPtr();

		if (atom(i)->getModel()->isBond())
		{
			bond = ToBondPtr(atom(i)->getModel());
		}
		else
		{
			continue;
		}

		Bond::setDampening(&*bond, -0.0);

		double kick = 0.6;
		std::string id = getMonomer()->getIdentifier();

		if (id == "tyr" || id == "phe" || id == "trp" || id == "his")
		{
			kick = 0.3;
		}

		if (bond->isRefinable() && atom(i)->getAtomName() == "CB")
		{
			Bond::setTorsionBlur(&*bond, kick);
		}
	}
}

void Sidechain::splitConformers(int count)
{
	if (count < 0)
	{
		count = conformerCount();
	}

	if (count <= 1) return;

	if (getMonomer()->getBackbone()->findAtoms("N").size() != 1)
	{
		std::cout << "Not splitting whole residue conformer, "
		<< getMonomer()->getResidueNum() << getMonomer()->getIdentifier() << std::endl;
//		return;
	}

	AtomPtr start = findAtom("CB");

	if (!start || !start->getModel()->isBond())
	{
		return;
	}

	BondPtr bond = ToBondPtr(start->getModel());

	for (int i = 1; i < count; i++)
	{
		bond->splitBond();
	}

	AtomList atoms = findAtoms("CB");

	for (int i = 0; i < atoms.size(); i++)
	{
		if (atoms[i].expired()) continue;

		AtomPtr atom = atoms[i].lock();
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

			if (atoms[i].expired())
			{
				continue;
			}

			AtomPtr atom = atoms[i].lock();
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
		if (cbs[i].expired()) continue;

		AtomPtr cb = cbs[i].lock();
		if (!cb->getModel()->isBond()) continue;

		BondPtr bond = ToBondPtr(cb->getModel());

		double occ = rotamers[i].allOccupancy;
		Bond::setOccupancy(&*bond, occ);
	}

	refreshRotamers();
}

void Sidechain::refreshRotamers()
{
	if (!_rotamerised)
	{
		warn_user("Trying to refresh rotamers without having rotamerised");
	}

	AtomList atoms = findAtoms("CB");
	double occTotal = 0;
	double occWeight = 0;
	double expTotal = 0;

	for (int i = 0; i < atoms.size(); i++)
	{
		AtomPtr atom = atoms[i].lock();
		if (atom->getModel()->isBond())
		{
			BondPtr bond = ToBondPtr(atom->getModel());
			double occ = Bond::getOccupancy(&*bond);
			double myExp = exp(-_exponent * occTotal * occTotal);

			expTotal += myExp;
			occTotal += occ;
			occWeight += occ * myExp;
		}
	}

	occWeight /= expTotal;

	occTotal = 0;
	double test = 0;
	for (int i = 0; i < atoms.size(); i++)
	{
		AtomPtr atom = atoms[i].lock();
		if (atom->getModel()->isBond())
		{
			BondPtr bond = ToBondPtr(atom->getModel());
			double myExp = exp(-_exponent * occTotal * occTotal);
			double mult = myExp / (occWeight * expTotal);
			occTotal += Bond::getOccupancy(&*bond);
			bond->setOccupancyMult(mult);

			double newOcc = bond->getMultOccupancy();
		//	std::cout << "New occupancy state " << i << ", " << newOcc << std::endl;
			test += newOcc;
		}
	}
}
