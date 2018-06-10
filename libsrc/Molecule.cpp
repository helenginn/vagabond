//
//  Molecule.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "shared_ptrs.h"
#include "Bucket.h"
#include "Plucker.h"
#include "Molecule.h"
#include "Atom.h"
#include "Absolute.h"
#include "Element.h"
#include "Bond.h"
#include <float.h>
#include <iostream>
#include <fstream>
#include "fftw3d.h"
#include "mat3x3.h"
#include "Options.h"
#include "CSV.h"
#include "Sidechain.h"
#include "Monomer.h"

Molecule::Molecule()
{
	_absoluteBFacSubtract = 0.0;
	_absoluteBFacMult = 1.0;
	_magicRotAxis = make_vec3(1, 0, 0);
	_rotationAxis = make_vec3(1, 0, 0);
	_rotationCentre = make_vec3(nan(" "), nan(" "), nan(" "));
	_sphereDiffOffset = empty_vec3();
	_rotationAngle = 0;
	_changedRotations = true;
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

void Molecule::tieAtomsUp()
{
	if (getClassName() != "Polymer")
	{
		double mult = Options::getBMult();
		std::cout << "Setting HETATM B factor multiplier to " << mult <<
		" for Chain " << getChainID() << std::endl;
		setAbsoluteBFacMult(mult);
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
	std::cout << "| Atoms: " << atomCount() << std::endl;
}

void Molecule::refine(CrystalPtr, RefinementType)
{

}

std::string Molecule::makePDB(PDBType, CrystalPtr)
{
	return "";
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
			some += atom(i)->getElement()->electronCount();
			someCount++;
		}

		total += atom(i)->getElement()->electronCount();
		totalCount++;
	}

	*tied += some;
	*all += total;
}

void Molecule::addAtom(AtomPtr atom)
{
	AtomGroup::addAtom(atom);
}


void Molecule::resetInitialPositions()
{
	for (int i = 0; i < atomCount(); i++)
	{
		vec3 pos = atom(i)->getAbsolutePosition();
		atom(i)->setInitialPosition(pos);

		ModelPtr model = atom(i)->getModel();

		if (model->isBond())
		{
			ToBondPtr(model)->resetBondDirection();
		}
	}
}

std::vector<mat3x3> Molecule::getExtraRotations()
{
	if (!_changedRotations && _extraRotationMats.size())
	{
		return _extraRotationMats;
	}

	_extraRotationMats.clear();

	calculateExtraRotations();

	_changedRotations = false;

	return _extraRotationMats;
}

void Molecule::addProperties()
{
	addStringProperty("chain_id", &_chainID);
	addDoubleProperty("absolute_bfac_mult", &_absoluteBFacMult);
	addDoubleProperty("absolute_bfac_subtract", &_absoluteBFacSubtract);
	addVec3ArrayProperty("centroids", &_centroids);
	addVec3ArrayProperty("centroid_offsets", &_centroidOffsets);
	addVec3ArrayProperty("trans_tensor_offsets", &_transTensorOffsets);
	addVec3Property("magic_rot_axis", &_magicRotAxis);
	addVec3Property("rotation_axis", &_rotationAxis);
	addVec3Property("rot_centre", &_rotationCentre);
	addDoubleProperty("rotation_angle", &_rotationAngle);
	addDoubleProperty("trans_exponent", &_transExponent);
	addMat3x3ArrayProperty("extra_rotations", &_extraRotationMats);
	addMat3x3ArrayProperty("rotations", &_rotations);
	
	exposeFunction("set_absolute_bfac_mult", vsSetAbsoluteBFacMult);
	exposeFunction("set_absolute_bfac_subtract", vsSetAbsoluteBFacSubtract);

	for (int i = 0; i < atomCount(); i++)
	{
		addChild("atom", atom(i));
	}
}

void Molecule::setAbsoluteBFacMult(double mult)
{
	_absoluteBFacMult = mult;
	std::cout << "Setting absolute B factor multiplier to " << mult << std::endl;
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


