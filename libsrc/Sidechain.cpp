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


void Sidechain::refine(CrystalPtr target, RefinementType rType)
{
	int resNum = getMonomer()->getResidueNum();
/*
	for (int i = 0; i < bondCount(); i++)
	{
		if (bond(i)->isNotJustForHydrogens() && bond(i)->isUsingTorsion())
		{
			if (rType == RefinementFine)
			{
				addTorsionBlur(bond(i), deg2rad(2.0), deg2rad(0.2));
			}
		}
		else if (bond(i)->isNotJustForHydrogens())
		{
			if (rType == RefinementFine)
			{
			}
		}
	}*/

	const char *atoms[] = {"CB", "CG", "CD", "CE"};
	std::vector<std::string> atomStrs(atoms, atoms+4);
/*
	if (resNum == 104)
	{
		AtomPtr myAtom = findAtom("CB");
		BondPtr bond = std::static_pointer_cast<Bond>(myAtom->getModel());

		setupGrid();

		addTorsion(bond, deg2rad(360.0), deg2rad(4.0));
		addTorsionBlur(bond, deg2rad(18.0), deg2rad(0.6));

		for (int j = 0; j < bond->downstreamAtomCount(); j++)
		{
			addSampled(bond->downstreamAtom(j));
		}

		setMock();
		setJobName("thr_trial_" + i_to_str(resNum));
		setCrystal(target);
		sample();
	}*/

	for (int i = 0; i < atomStrs.size(); i++)
	{
		std::string atom = atomStrs[i];
		AtomPtr myAtom = findAtom(atom);

		if (!myAtom)
		{
			continue;
		}

		BondPtr bond = std::static_pointer_cast<Bond>(myAtom->getModel());

		if (!bond->isUsingTorsion())
		{
			continue;
		}

		std::string preAtom = bond->getMajor()->getAtomName();

		if (preAtom != "CA")
		{
			setupNelderMead();
			addTorsionNextBlur(bond, 0.1, deg2rad(0.1));

			for (int j = 0; j < bond->downstreamAtomCount(); j++)
			{
				addSampled(bond->downstreamAtom(j));
			}

			setJobName("compensate_" + preAtom + "_" + atom + "_" + i_to_str(resNum));
			setCrystal(target);
			sample();
		}
/*
		if (preAtom == "CB")
		{
			ModelPtr preModel = bond->getParentModel();
			if (preModel->getClassName() != "Bond")
			{
				break;
			}

			BondPtr preBond = std::static_pointer_cast<Bond>(preModel);

			setupGrid();
			addTorsion(preBond, deg2rad(180.0), deg2rad(4.0));
			addBendAngle(bond, deg2rad(90.0), deg2rad(1.0));

			addSampled(bond->getMinor());

			setJobName("grid_" + preAtom + "_" + atom + "_" + i_to_str(resNum));
			setCrystal(target);
			sample();
		}
*/
		ModelPtr preModel = bond->getParentModel();

		if (preModel->getClassName() == "Bond")
		{
			setupNelderMead();

			BondPtr preBond = std::static_pointer_cast<Bond>(preModel);

			addBendAngle(bond, deg2rad(0.5), deg2rad(0.1));
			addBendBlur(bond, deg2rad(0.3), deg2rad(0.1));

			addSampled(bond->getMinor());

			setJobName("bend_" + preAtom + "_" + atom + "_" + i_to_str(resNum));
			setCrystal(target);
			sample();
		}

		setupNelderMead();

		addTorsion(bond, deg2rad(1.0), deg2rad(0.2));
		addTorsionBlur(bond, deg2rad(1.5), deg2rad(0.2));

		for (int j = 0; j < bond->downstreamAtomCount(); j++)
		{
			addSampled(bond->downstreamAtom(j));
		}

		setJobName("torsion_" + preAtom + "_" + atom + "_" + i_to_str(resNum));
		setCrystal(target);
		sample();


	}


/*
	if (!sampleSize())
	{
		for (int i = 0; i < atomCount(); i++)
		{
			if (atom(i)->getAtomName() != "CA")
			{
				addSampled(atom(i));
			}
		}
	}
*/
//	setCrystal(target);
//	sample();

}