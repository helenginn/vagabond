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

#include <hcsrc/FileReader.h>
#include "TotalModel.h"
#include "SpaceSample.h"
#include "WaterNetwork.h"
#include "ConfSpace.h"
#include "Options.h"
#include "Shouter.h"
#include "Polymer.h"
#include "Motion.h"
#include "Anchor.h"
#include "Atom.h"

TotalModel::TotalModel()
{
	_tied = false;
	_sampleNum = -1;
	_confSpace = NULL;
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

void TotalModel::refreshAnchors(bool space)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->isPolymer())
		{
			PolymerPtr pol = ToPolymerPtr(molecule(i));

			if (pol->getAnchorModel())
			{
				pol->getAnchorModel()->forceRefresh();
				SpaceSample *sp = pol->getAnchorModel()->spaceSample();
				if (sp && space)
				{
					sp->generatePoints(shared_from_this(), true);
				}
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
	addIntProperty("sample_num", &_sampleNum);
	addBoolProperty("tied", &_tied);
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

MotionPtr TotalModel::getOverallMotion()
{
	if (_motions.size())
	{
		for (int i = 0; i < motionCount(); i++)
		{
			if (_motions[i]->getName() == "all")
			{
				return _motions[i];
			}
		}
	}

	return MotionPtr();
}

void TotalModel::makeOverallMotion()
{
	if (!_tied)
	{
		return;
	}

	if (getOverallMotion() != MotionPtr())
	{
		return;
	}
	
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer() ||
		    !ToPolymerPtr(molecule(i))->getAnchorModel())
		{
			continue;
		}

		PolymerPtr pol = ToPolymerPtr(molecule(i));
		pol->getAnchorModel()->atLeastOneMotion();
	}
}

int TotalModel::getSampleNum()
{
	if (Options::getNSamples() >= 0)
	{
		_sampleNum = Options::getNSamples();
		Options::setNSamples(NULL, -1);
	}

	if (_sampleNum < 0) 
	{
		_sampleNum = 120;
	}

	double totalPoints = _sampleNum;
	
	if (totalPoints < 0)
	{
		totalPoints = 0;
	}
	
	return totalPoints;
}

void TotalModel::tiedUpScattering()
{
	double tied = 0;
	double total = 0;
	int flex = 0;
	int pos = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tiedUpScattering(&tied, &total);
		molecule(i)->reportParameters();
		molecule(i)->addParamCounts(&pos, &flex);
	}
	
	std::cout << "Total positional params: " << pos << std::endl;
	std::cout << "Total flexibility params: " << flex << std::endl;
	std::cout << "Total params: " << pos + flex << std::endl;

	std::cout << std::fixed << std::setprecision(0);
	std::cout << "Tied up " << 100. * sqrt(tied / total) << "% of"\
	" the scattering electrons." << std::endl;
	std::cout << std::endl;
}

double TotalModel::averageBFactor()
{
	double ave = 0;
	double count = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		MoleculePtr mole = molecule(i);

		if (!mole->isPolymer())
		{
			continue;
		}
		
		ave += ToPolymerPtr(mole)->getAverageBFactor();
		count++;
	}
	
	return ave/count;
}

WaterNetworkPtr TotalModel::getWaterNetwork()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isWaterNetwork())
		{
			continue;
		}
		
		WaterNetworkPtr wat = ToWaterNetworkPtr(molecule(i));
		return wat;
	}

	return WaterNetworkPtr();
}

void TotalModel::resetMotions()
{
	for (int i = 0; i < _motions.size(); i++)
	{
		_motions[i]->reset();
	}
	
	refreshPositions();
}

void TotalModel::setupConformationalSpace()
{
	if (_confSpace != NULL)
	{
		return;
	}
	
	std::string filename = "membership.csv";

	_confSpace = new ConfSpace(16);

	for (size_t i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		PolymerPtr pol = ToPolymerPtr(molecule(i));
		ConfSpace *confSpace = new ConfSpace(16);
		_spaces[pol] = confSpace;

		if (file_exists(filename))
		{
			confSpace->readFromFile(filename);
		}
		else
		{
			confSpace->calculateFrom(pol);
		}
		
		AnchorPtr anchor = pol->getAnchorModel();

		anchor->makeSpaceSample(confSpace);
		SpaceSample *sp = anchor->spaceSample();
		AtomList cas = pol->findAtoms("CA");
		sp->setAtoms(cas);
	}
	
	for (size_t i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}
		
		PolymerPtr pol = ToPolymerPtr(molecule(i));
		AnchorPtr anchor = pol->getAnchorModel();
	}
}
