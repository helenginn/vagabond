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

bool Sidechain::shouldRefineMagicAxis(BondPtr bond)
{
	if (_refinedMagicAxisCount > 2) return false;

	return (bond->getMinor()->getAtomName() == "CB");
}


void Sidechain::refine(CrystalPtr target, RefinementType rType)
{
	int resNum = getMonomer()->getResidueNum();

	if (!isTied())
	{
		return;
	}

	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr myAtom = atom(i);
		ModelPtr model = myAtom->getModel();

		if (!model->isBond())
		{
			continue;
		}

		BondPtr bond = std::static_pointer_cast<Bond>(model);

		if (!bond->isRefinable())
		{
			continue;
		}

		int groups = bond->downstreamAtomGroupCount();

		for (int k = 0; k < groups; k++)
		{
			if (rType != RefinementModelRMSD)
			{
				break;
			}

			if (shouldRefineMagicAxis(bond))
			{
				_refinedMagicAxisCount++;

				/* Find the vague direction */
				setupGrid();
				addSampledAtoms(shared_from_this());
				addMagicAxisBroad(bond);
				setSilent();
				setJobName("broad_axis_" +  bond->shortDesc());
				setScoreType(ScoreTypeModelRMSDZero);
				sample();

				/* Fine-tune the axis */
				setupNelderMead();
				addMagicAxis(bond, deg2rad(10.0), deg2rad(2.0));
				setJobName("magic_axis_" +  bond->shortDesc());
				setSilent();
				addSampledAtoms(shared_from_this());
				setScoreType(ScoreTypeModelRMSDZero);
				sample();

				bond->resetAxis();
			}

			/* Refine model position */

			setupNelderMead();
			setupTorsionSet(bond, k, 3, resNum, ANGLE_SAMPLING, deg2rad(0.01));

			setScoreType(ScoreTypeModelPos);
			setSilent();
			setJobName("model_pos_" +  bond->shortDesc());
			sample();
		}

		if (rType != RefinementFineBlur)
		{
			continue;
		}

		for (int k = 0; k < groups; k++)
		{
			setupNelderMead();
			setupTorsionSet(bond, k, 5, resNum, 0.05, 0.2);
			setCrystal(target);
			sample();

			setupNelderMead();
			setupTorsionSet(bond, k, 4, resNum, 0.3, 0.02);
			addDampening(bond, 0.02, 0.1);
			setCrystal(target);
			sample();
		}


		bond->setActiveGroup(0);

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

