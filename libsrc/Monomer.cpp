//
//  Monomer.cpp
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Monomer.h"
#include "Backbone.h"
#include "Sidechain.h"
#include "Atom.h"
#include "Knotter.h"
#include "Bond.h"
#include "shared_ptrs.h"
#include "Polymer.h"
#include "Shouter.h"
#include "Options.h"

Monomer::Monomer()
{
    _backbone = BackbonePtr(new Backbone());
    _sidechain = SidechainPtr(new Sidechain());
}

void Monomer::setup()
{
    _backbone->setMonomer(shared_from_this());
    _backbone->setPolymer(getPolymer());
    _backbone->setResNum(_residueNum);
    _sidechain->setMonomer(shared_from_this());
    _sidechain->setResNum(_residueNum);

}

void Monomer::addAtom(AtomPtr atom)
{
    std::string newAtomName = atom->getAtomName();

    AtomPtr existingAtom = findAtom(newAtomName);

    if (existingAtom)
    {
    //    std::string atomDesc = existingAtom->shortDesc();
    //    warn_user("There is an existing conformation for " + atomDesc
    //              + "\nCurrently only side chains are supported.");
    }

    AtomGroup::addAtom(atom);

    atom->setMonomer(shared_from_this());
    
    bool isBoth = atom->isBackboneAndSidechain();

    if (isBoth)
    {
        if (!existingAtom)
        {
            _backbone->addAtom(atom);
        }

        _sidechain->addAtom(atom);
        return;
    }

    bool isBackbone = atom->isBackbone();

    if (isBackbone && !existingAtom)
    {
        _backbone->addAtom(atom);
    }
    else if (!isBackbone)
    {
        _sidechain->addAtom(atom);
    }
}

void Monomer::setSideKick(double value)
{
    AtomPtr atom = findAtom("CB");

    if (!atom)
    {
        return;
    }

    ModelPtr model = atom->getModel();
    if (!model->isBond())
    {
        return;
    }

    BondPtr bond = ToBondPtr(model);

    if (bond->isFixed() || !bond->isUsingTorsion())
    {
        return;
    }

    Bond::setTorsionBlur(&*bond, value);
}

void Monomer::setKick(double value, bool beforeAnchor)
{
    AtomPtr atom;

    if (beforeAnchor)
    {
        atom = findAtom("C");
    }
    else
    {
        atom = findAtom("CA");
    }

    ModelPtr model = atom->getModel();
    if (model->getClassName() != "Bond")
    {
        shout_at_helen("Can't kick something that\n isn't a bond!");
    }

    BondPtr bond = ToBondPtr(model);

    Bond::setTorsionBlur(&*bond, value);
}

double Monomer::getKick()
{
    AtomPtr ca = findAtom("CA");
    ModelPtr model = ca->getModel();
    if (model->getClassName() != "Bond")
    {
        shout_at_helen("Can't kick something that\n isn't a bond!");
    }

    BondPtr bond = ToBondPtr(model);
    
    return Bond::getTorsionBlur(&*bond);
}

void Monomer::tieAtomsUp()
{
    KnotterPtr knotter = KnotterPtr(new Knotter());

    const int start = getPolymer()->getAnchor();

    knotter->setBackbone(_backbone);

    if (getResidueNum() >= start)
    {
        knotter->tieTowardsCTerminus();
    }
    else
    {
        knotter->tieTowardsNTerminus();
    }

    knotter->setSidechain(_sidechain);
    knotter->tie();

    if (getResidueNum() == start)
    {
        BondPtr bond = ToBondPtr(getBackbone()->findAtom("CA")->getModel());
    }
    else if (getResidueNum() == start - 1)
    {
        BondPtr bond = ToBondPtr(getBackbone()->findAtom("C")->getModel());
    }

    _backbone->setTied();
    _sidechain->setTied();
}

void Monomer::setBackboneDampening(double value)
{
    for (int i = 0; i < atomCount(); i++)
    {
        if (atom(i)->getModel()->isBond())
        {
            BondPtr bond = ToBondPtr(atom(i)->getModel());

            if (bond->isRefinable())
            {
                Bond::setDampening(&*bond, value);
            }
        }
    }
}

void Monomer::setSidechainDampening(double value)
{
    for (int i = 0; i < getSidechain()->atomCount(); i++)
    {
        ModelPtr model = getSidechain()->atom(i)->getModel();

        if (model->isBond())
        {
            Bond::setDampening(&*model, value);
        }
    }
}

bool Monomer::isAfterAnchor()
{
    int anchor = getPolymer()->getAnchor();
    bool after = (_residueNum >= anchor);

    return after;
}

std::string Monomer::getResCode()
{
    return GeomTable::getResCode(getIdentifier());
}

