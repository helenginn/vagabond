//  Crystal.h
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//
//

#include "shared_ptrs.h"
#include "Bucket.h"
#include "Molecule.h"
#include "Atom.h"
#include "Absolute.h"
#include "Bond.h"
#include <float.h>
#include <iostream>
#include <fstream>
#include "mat3x3.h"
#include "Options.h"
#include "CSV.h"
#include "Chelate.h"
#include "Sidechain.h"
#include "Monomer.h"
#include "Element.h"

Molecule::Molecule()
{
	_absoluteBFacSubtract = 0.0;
	_absoluteBFacMult = 1.0;
}

double Molecule::getAbsoluteBFacMult()
{
	_absoluteBFacMult = Options::getBMult();
	return _absoluteBFacMult;
}

double Molecule::getAbsoluteBFacSubt()
{
	if (Options::getBSubt() > 0)
	{
		_absoluteBFacSubtract = Options::getBSubt();
	}

	return _absoluteBFacSubtract;
}

void Molecule::tieAtomsUp()
{
	if (getClassName() != "Polymer")
	{
		double mult = Options::getBMult();
		setAbsoluteBFacMult(mult);
	}
}

void Molecule::summary()
{
	std::cout << "| I am chain " << getChainID() << std::endl;
	
	int hydrogen = 0;
	for (int i = 0; i < atomCount(); i++)
	{
		if (atom(i)->getElectronCount() == 1)
		{
			hydrogen++;
		}
	}
	
	std::cout << "| Atoms: " << atomCount() - hydrogen;
	
	if (hydrogen > 0)
	{
		std::cout << " (+ " << hydrogen << " hydrogens)";
	}
	
	std::cout << std::endl;
}

std::string Molecule::makePDB(PDBType pdbType, CrystalPtr crystal, 
                              int conformer)
{
	std::ostringstream stream;
	stream << getPDBContribution(pdbType, crystal, conformer);

	return stream.str();
}

void Molecule::reportParameters()
{

}

void Molecule::tiedUpScattering(double *tied, double *all)
{
	double total = 0;
	double totalCount = 0;
	double some = 0;
	double someCount = 0;

	for (int i = 0; i < atomCount(); i++)
	{
		if (!atom(i)->getElement())
		{
			std::cout << "Warning! Atom has no element: " << atom(i)->shortDesc() << std::endl;
		}

		if (!atom(i)->getModel())
		{
			std::cout << "Warning! Atom has no model: " << atom(i)->shortDesc() << std::endl;
		}

		if (atom(i)->getModel()->isBond())
		{
			some += atom(i)->getElectronCount();
			someCount++;
		}

		total += atom(i)->getElectronCount();
		totalCount++;
	}

	*tied += some;
	*all += total;
}

void Molecule::addAtom(AtomPtr atom)
{
	AtomGroup::addAtom(atom);
	Options::getActiveCrystal()->addAtom(atom);
}

void Molecule::addProperties()
{
	addStringProperty("chain_id", &_chainID);
	addDoubleProperty("bfac_subt", &_absoluteBFacSubtract);
	addDoubleProperty("bfac_mult", &_absoluteBFacMult);
	
	for (int i = 0; i < atomCount(); i++)
	{
		addChild("atom", atom(i));
	}
}

void Molecule::setAbsoluteBFacMult(double mult)
{
	_absoluteBFacMult = mult;
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	crystal->addComment("Absolute B factor multiplier changed to "
	+ f_to_str(mult, 2) + " for " + getChainID());
}

void Molecule::postParseTidy()
{
	for (int i = 0; i < atomCount(); i++)
	{
		ModelPtr model = atom(i)->getModel();
		if (!model) 
		{
			continue;
		}

		model->setMolecule(shared_from_this());
	}
}

std::vector<AtomPtr> Molecule::getCloseAtoms(AtomPtr one, double tol, bool cache)
{
	std::vector<AtomPtr> atoms;
	std::vector<AtomPtr> *searchAtoms = &_atoms;
	std::vector<AtomPtr> *atomListPtr = &atoms;
	
	if (cache)
	{
		atomListPtr = &_closeishAtoms;	
	}
	
	if (!cache && _closeishAtoms.size())
	{
		searchAtoms = &_closeishAtoms;
	}

	for (int i = 0; i < searchAtoms->size(); i++)
	{
		AtomPtr search = searchAtoms->at(i);
		if (one == search)
		{
			continue;
		}
		
		if (search->getMonomer())
		{
			SidechainPtr side = search->getMonomer()->getSidechain();
			if (side->isRotamerised())
			{
				continue;	
			}
		}

		if (!one->closeToAtom(search, tol))
		{
			continue;
		}
		
		std::vector<AtomPtr>::iterator it;
		it = std::find(atomListPtr->begin(), atomListPtr->end(), search);

		if (it != atomListPtr->end())
		{
			continue;	
		}
		
		atomListPtr->push_back(search);
	}
	
	return atoms;
}


void Molecule::setAbsoluteBFacSubtract(void *object, double subtract)
{
	Molecule *obj = static_cast<Molecule *>(object);
	obj->_absoluteBFacSubtract = subtract;
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	crystal->addComment("Absolute B factor subtractor changed to "
	+ f_to_str(subtract, 2) + " for " + obj->getChainID());
}

void Molecule::chelate(std::string element, double bufferB)
{
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();

	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr a = atom(i);
		
		if (a->getElement()->getSymbol() != element)
		{
			continue;
		}
		
		std::vector<AtomPtr> atoms;
		atoms = crystal->getCloseAtoms(a, 2.6);
		ChelatePtr chelate = ChelatePtr(new Chelate());
		chelate->setChelatedAtom(a);
		
		if (atoms.size() == 0)
		{
			continue;
		}
		
		for (size_t j = 0; j < atoms.size(); j++)
		{
			AtomPtr b = atoms[j];
			
			if (b == a)
			{
				continue;
			}
			
			chelate->addChelatingAtom(b);
		}
		
		chelate->setBufferB(bufferB);
		
		std::cout << "Chelated " << chelate->shortDesc() << std::endl;
	}
}

void Molecule::refitToSavedPositions()
{
	std::cout << "Dud." << std::endl;
}
