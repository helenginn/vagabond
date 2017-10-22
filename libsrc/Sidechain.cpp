//
//  Sidechain.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Sidechain.h"
#include "Sampler.h"
#include "Bond.h"
#include "Atom.h"
#include <iostream>
#include "Monomer.h"
#include "FileReader.h"
#include "Absolute.h"

void Sidechain::refine(CrystalPtr target, RefinementType rType)
{
	int resNum = getMonomer()->getResidueNum();

	if (!isTied())
	{
		return;
	}

	_timesRefined++;

	const char *atoms[] = {"CB", "OG", "CG", "SD", "CG1", "CG2", "CD", "CE", "NZ", "OG1"};
	std::vector<std::string> atomStrs(atoms, atoms+8);

	for (int i = 0; i < atomStrs.size(); i++)
	{
		std::string atom = atomStrs[i];
		AtomList myAtoms = findAtoms(atom);

		if (!myAtoms.size())
		{
			continue;
		}

		for (int j = 0; j < myAtoms.size(); j++)
		{
			AtomPtr myAtom = myAtoms[j].lock();
			ModelPtr model = myAtom->getModel();

			BondPtr bond = std::static_pointer_cast<Bond>(model);
			std::string majorAtom = bond->getMajor()->getAtomName();

			if (!bond->isNotJustForHydrogens() || bond->isFixed())
			{
				continue;
			}

			ModelPtr preModel = bond->getParentModel();

			int groups = bond->downstreamAtomGroupCount();

			for (int k = 0; k < groups; k++)
			{
				if (rType != RefinementModelRMSD)
				{
					break;
				}

				ScoreType scoreType = ScoreTypeModelRMSDZero;
				for (int l = 0; l < 1; l++)
				{
					if (bond->getMinor()->getAtomName() != "CB")
					{
						continue;
					}

					if (scoreType == ScoreTypeModelRMSDZero)
					{
						setupGrid();
						addSampledAtoms(shared_from_this());
						addMagicAxisBroad(bond);
						setSilent();
						setJobName("broad_axis_" +  bond->shortDesc());
						setScoreType(scoreType);
						sample();
					}

					setupNelderMead();

					if (scoreType == ScoreTypeModelRMSD)
					{
						addDampening(bond, 0.05, 0.05);
					}
					else
					{
						addMagicAxis(bond, deg2rad(10.0), deg2rad(2.0));
					}

					setJobName("magic_axis_" +  bond->shortDesc());
					setSilent();
					addSampledAtoms(shared_from_this());
					setScoreType(scoreType);
					sample();

					if (scoreType == ScoreTypeModelRMSD)
					{
						scoreType = ScoreTypeModelRMSDZero;
					}
					else
					{
						scoreType = ScoreTypeModelRMSD;

					}

					bond->resetAxis();
				}
/*
				rType = RefinementModelPos;
			}

			if (rType == RefinementModelPos)
			{*/
				for (int k = 0; k < groups; k++)
				{
					setupNelderMead();
					setupTorsionSet(bond, k, 3, resNum, ANGLE_SAMPLING, deg2rad(0.01));
					setScoreType(ScoreTypeModelPos);
					setSilent();
					setJobName("model_pos_" +  bond->shortDesc());
					sample();
					rType = RefinementModelPos;
				}

			}

			if (rType != RefinementFineBlur)
			{
				return;
			}
			
			for (int k = 0; k < groups; k++)
			{
				setupNelderMead();
				setupTorsionSet(bond, k, 5, resNum, 2.5, 0.1);
				setCrystal(target);
				sample();

				setupNelderMead();
				setupTorsionSet(bond, k, 5, resNum, 0.5, 0.1);
				addDampening(bond, 0.02, 0.1);
				setCrystal(target);
				sample();
			}


			bond->setActiveGroup(0);

		}

	}
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

