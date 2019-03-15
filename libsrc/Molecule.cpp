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
#include "fftw3d.h"
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

void Molecule::makePowderList()
{
	CSVPtr csvAngles = CSVPtr(new CSV(3, "dist1", "dist2", "angle"));
	CSVPtr csvDist = CSVPtr(new CSV(1, "dist1"));
	double maxDistance = 5.0;

	for (int i = 1; i < atomCount() - 1; i++)
	{
		vec3 iPos = atom(i)->getAbsolutePosition();

		for (int j = 0; j < i; j++)
		{
			vec3 jPos = atom(j)->getAbsolutePosition();

			vec3 diff1 = vec3_subtract_vec3(jPos, iPos);

			double aLength = vec3_length(diff1);
			if (aLength > maxDistance)
			{
				continue;
			}

			csvDist->addEntry(1, aLength);

			for (int k = 0; k < atomCount(); k++)
			{
				if (k == j) continue;
				if (k == i) continue;

				vec3 kPos = atom(k)->getAbsolutePosition();
				vec3 diff2 = vec3_subtract_vec3(kPos, jPos);
				double cosine = vec3_cosine_with_vec3(diff2, diff1);
				double angle = acos(cosine); 
				//                if (cosine < 0) angle += deg2rad(90);

				if (angle != angle) continue;

				double bLength = vec3_length(diff2);

				if (bLength > maxDistance)
				{
					continue;
				}

				csvAngles->addEntry(3, aLength, bLength, rad2deg(angle));
				csvAngles->addEntry(3, bLength, aLength, rad2deg(angle));
				vec3 diff3 = vec3_subtract_vec3(kPos, iPos);
				cosine = vec3_cosine_with_vec3(diff3, diff1);
				angle = acos(cosine); 
				//                if (cosine < 0) angle += deg2rad(90);

				if (angle != angle) continue;

				double cLength = vec3_length(diff3);
				if (cLength > maxDistance)
				{
					continue;
				}

				csvAngles->addEntry(3, aLength, cLength, rad2deg(angle));
				csvAngles->addEntry(3, cLength, aLength, rad2deg(angle));
			}
		}
	}

	csvDist->writeToFile(getChainID() + "_distances.csv");
	csvAngles->writeToFile(getChainID() + "_angles.csv");
}

void Molecule::forceModelRecalculation()
{
	for (int i = 0; i < atomCount(); i++)
	{
		ModelPtr model = atom(i)->getModel();
		model->recalculate();
	}
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

void Molecule::addToSolventMask(FFTPtr fft, mat3x3 _real2frac, double rad,
	                           std::vector<Atom *> *ptrs, int conf)
{
	vec3 offset = make_vec3(0, 0, 0);

	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->addToSolventMask(fft, _real2frac, rad, ptrs, conf);
	}
}

void Molecule::addToMap(FFTPtr fft, mat3x3 _real2frac, bool mask)
{
	vec3 offset = make_vec3(0, 0, 0);

	for (int i = 0; i < atomCount(); i++)
	{
		atom(i)->addToMap(fft, _real2frac, offset, mask);
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

void Molecule::refine(CrystalPtr crystal, RefinementType type)
{
	for (int i = 0; i < atomCount(); i++)
	{
		AtomPtr anAtom = atom(i);
		if (anAtom->getElement()->getSymbol() != "O" || 
		    !anAtom->isHeteroAtom())
		{
			continue;
		}
		
		ModelPtr model = anAtom->getModel();
		if (!model || !model->isAbsolute())
		{
			continue;
		}

		AbsolutePtr abs = ToAbsolutePtr(model);
		
		setupStepSearch();
		setCrystal(crystal);
		setCycles(20);
		addSampled(anAtom);	
//		addAbsolutePosition(abs, 0.005, 0.0001);
		setJobName("xyz_" + anAtom->shortDesc());
		sample();
	}
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

	refreshBModels();
}

void Molecule::refreshBModels()
{
	for (int i = 0; i < atomCount(); i++)
	{
		ModelPtr model = atom(i)->getModel();
		model->recalculate();
	}
	propagateChange();
	refreshPositions();
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
	obj->refreshBModels();
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

void Molecule::rigidBodyFit()
{
	std::cout << "Dud." << std::endl;
}
