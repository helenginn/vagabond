//
//  Molecule.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "shared_ptrs.h"
#include "Molecule.h"
#include "Atom.h"
#include "Element.h"
#include "Bond.h"
#include <float.h>
#include <iostream>
#include "fftw3d.h"
#include "mat3x3.h"
#include "Options.h"

Molecule::Molecule()
{
    _absoluteBFacSubtract = 0.0;
    _absoluteBFacMult = 1.0;
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

void Molecule::addToMap(FFTPtr fft, mat3x3 _real2frac)
{
    for (int i = 0; i < atomCount(); i++)
    {
        atom(i)->addToMap(fft, _real2frac);
    }
}

void Molecule::summary()
{
    std::cout << "| I am chain " << getChainID() << std::endl;
    std::cout << "| Atoms: " << atomCount() << std::endl;
}

void Molecule::refine(CrystalPtr target, RefinementType rType)
{

}

std::string Molecule::makePDB(PDBType pdbType, CrystalPtr crystal)
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
        if (atom(i)->getModel()->getClassName() == "Bond")
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

void Molecule::resetInitialPositions()
{
    for (int i = 0; i < atomCount(); i++)
    {
        vec3 pos = atom(i)->getPosition();
        atom(i)->setInitialPosition(pos);

        ModelPtr model = atom(i)->getModel();

        if (model->isBond())
        {
            ToBondPtr(model)->resetBondDirection();
        }
    }
}


void Molecule::addProperties()
{
    addStringProperty("chain_id", &_chainID);
    addDoubleProperty("absolute_bfac_mult", &_absoluteBFacMult);
    addDoubleProperty("absolute_bfac_subtract", &_absoluteBFacSubtract);

    for (int i = 0; i < atomCount(); i++)
    {
        addChild("atom", atom(i));
    }
}

