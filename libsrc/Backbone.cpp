//
//  Backbone.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Backbone.h"
#include "Monomer.h"
#include "Bond.h"
#include "Atom.h"
#include "Absolute.h"
#include "FileReader.h"
#include "Sidechain.h"

bool Backbone::shouldRefineAtom(AtomPtr atom)
{
	return true;
}

void Backbone::refine(CrystalPtr target, RefinementType rType)
{
	if (!isTied())
	{
		return;
	}
	
	if (rType == RefinementSidechain)
	{
		return;
	}

	std::cout << getMonomer()->getResCode() << std::flush;

	if (!paramCount())
	{
		double range = 2.;

		if (_timesRefined >= 3)
		{
			range = 0.2;
		}
		
		if (rType == RefinementCrude)
		{
			range /= 2;
		}

		switch (rType)
		{
			case RefinementModelPos:
			addParamType(ParamOptionTorsion, range);
			addParamType(ParamOptionBondAngle, range / 1);
			addParamType(ParamOptionNumBonds, 3);
			break;
			
			case RefinementCrude:
			addParamType(ParamOptionTorsion, range);
			addParamType(ParamOptionTwist, range);
			addParamType(ParamOptionNumBonds, 5);
			break;
			
			case RefinementFine:
			addParamType(ParamOptionTorsion, range);
			addParamType(ParamOptionNumBonds, 4);
			break;

			case RefinementRMSDZero:
			addParamType(ParamOptionMagicAngles, 3);
			addParamType(ParamOptionNumBonds, 4);
			break;

			default:
			break;
		}
		
		if (_timesRefined >= 3)
		{
			switch (rType)
			{
				case RefinementModelPos:
				addParamType(ParamOptionBondAngle, range / 1);
				break;

				case RefinementFine:
				addParamType(ParamOptionBondAngle, 0.1);

				default:
				break;
			}
		}
	}

	MonomerPtr monomer = getMonomer();
	_includeForRefine = monomer->includingInRefinement();
	
	SidechainPtr sidechain = monomer->getSidechain();
	if (sidechain && rType == RefinementFine)
	{
//		addIncludeForRefinement(sidechain);
	}

	AtomGroup::refine(target, rType);
	clearParams();
	clearIncludeForRefinements();
}

AtomPtr Backbone::betaCarbonTorsionAtom()
{
	/* What is the major atom of the bond describing CA? */

	AtomPtr ca = findAtom("CA");

	if (!ca)
	{
		return AtomPtr();
	}

	ModelPtr model = ca->getModel();

	if (model->getClassName() == "Bond")
	{
		BondPtr bond = ToBondPtr(model);
		return bond->getMajor();
	}

	return AtomPtr();
}

void Backbone::addProperties()
{
	addIntProperty("res_num", &_resNum);
	AtomGroup::addProperties();
}
