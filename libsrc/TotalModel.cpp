// vagabond
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

#include "TotalModel.h"
#include "Options.h"
#include "Shouter.h"
#include "Polymer.h"
#include "Motion.h"
#include "Anchor.h"
#include "Atom.h"

TotalModel::TotalModel()
{

}

void TotalModel::removeMolecule(MoleculePtr mol)
{
	std::cout << "Removing molecule " << mol->getChainID() << std::endl;

	for (int i = 0; i < mol->atomCount(); i++)
	{
		if (mol->isPolymer())
		{
			PolymerPtr pol = ToPolymerPtr(mol);
			pol->removeAtom(mol->atom(i));
		}

		removeAtom(mol->atom(i));
		i--;
	}
	
	_molecules.erase(mol->getChainID());
}

void TotalModel::removeAtom(AtomPtr atom)
{
	std::cout << "Removing atom " << atom->longDesc() << std::endl;

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->removeAtom(atom);
	}
	
	AtomGroup::removeAtom(atom);
}

void TotalModel::refreshAnchors()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->isPolymer())
		{
			PolymerPtr pol = ToPolymerPtr(molecule(i));

			if (pol->getAnchorModel())
			{
				pol->getAnchorModel()->forceRefresh();
			}
		}
	}

	refreshPositions();
}

void TotalModel::setAnchors()
{
	std::string anchor = Options::anchorString();
	std::vector<std::string> components = split(anchor, ',');
	
	for (int i = 0; i < components.size(); i++)
	{
		std::string chain = "";
		std::string number = "";

		for (int j = 0; j < components[i].length(); j++)
		{
			if (components[i][j] >= '0' && components[i][j] <= '9')
			{
				number.push_back(components[i][j]);
			}
			else
			{
				chain.push_back(components[i][j]);
			}
		}
		
		int anchorPoint = atoi(number.c_str());
		
		for (int i = 0; i < moleculeCount(); i++)
		{
			std::string thisChain = molecule(i)->getChainID();
			std::string truncated = thisChain.substr(0, chain.length());
			
			if (molecule(i)->isPolymer() && 
			    chain == truncated)
			{
				PolymerPtr polymer = ToPolymerPtr(molecule(i));
				polymer->setAnchor(anchorPoint);
				std::cout << "Setting custom anchor " << anchorPoint;
				std::cout << " for chain " << chain << std::endl;
			}
		}

	}

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() == "Polymer")
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			if (polymer->getAnchor() < 0)
			{
				polymer->findAnchorNearestCentroid();
			}
		}
	}
}

size_t TotalModel::polymerCount()
{
	size_t count = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		count++;
	}

	return count;
}

void TotalModel::recalculateAtoms()
{
	empty();
	
	for (int i = 0; i < moleculeCount(); i++)
	{
		addAtomsFrom(molecule(i));
	}
}

void TotalModel::postParseTidy()
{
	for (int i = 0; i < motionCount(); i++)
	{
		_motions[i]->updateAtoms();
		_motions[i]->absorbScale();
	}

	recalculateAtoms();
}

void TotalModel::addProperties()
{
	addStringProperty("filename", &_filename);

	for (int i = 0; i < moleculeCount(); i++)
	{
		addChild("molecule", molecule(i));
	}
}

void TotalModel::addObject(ParserPtr object, std::string category)
{
	if (category == "molecule")
	{
		MoleculePtr molecule = ToMoleculePtr(object);
		addMolecule(molecule);
	}

	if (category == "motion")
	{
		MotionPtr motion = ToMotionPtr(object);
		addMotion(motion);
	}

}

void TotalModel::addMotion(MotionPtr mot, PolymerPtr origPol)
{
	_motions.push_back(mot);
	
	if (!origPol)
	{
		return;
	}
	
	char first = origPol->getChainID()[0];

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer() || molecule(i) == origPol)
		{
			continue;
		}
		
		if (!ToPolymerPtr(molecule(i))->isFullyTied())
		{
			continue;
		}
		
		if (molecule(i)->getChainID()[0] == first)
		{
			mot->addToPolymer(ToPolymerPtr(molecule(i)));
		}
	}
}

void TotalModel::addMolecule(MoleculePtr molecule)
{
	if (!molecule)
	{
		return;
	}

	if (molecule->getChainID().length() <= 0)
	{
		shout_at_helen("Polymer chain ID is missing while trying\n"\
		               "to interpret PDB file.");
	}

	_molecules[molecule->getChainID()] = molecule;
}
