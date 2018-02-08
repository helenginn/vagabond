//
//  Sidechain.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Sidechain__
#define __vagabond__Sidechain__

#include <stdio.h>
#include "AtomGroup.h"
#include "Sampler.h"
#include "FileReader.h"

class Sidechain : public AtomGroup
{
public:
    Sidechain()
    {
        _canRefine = false;
        _rotamerised = false;
        _exponent = 0;
    }

    bool canRefine()
    {
        return _canRefine;
    }

    void setCanRefine(bool canRefine)
    {
        _canRefine = canRefine;
    }

    void setResNum(int resNum)
    {
        _resNum = resNum;
    }

    void setPolymer(PolymerPtr poly)
    {
        _myPolymer = poly;
    }

    PolymerPtr getPolymer()
    {
        return _myPolymer.lock();
    }

    virtual bool shouldRefineAngles()
    {
        return (_timesRefined > 0);
    }

    void setInitialDampening();
    void fixBackboneTorsions(AtomPtr betaTorsion);
    void splitConformers(int count = -1);
    void parameteriseAsRotamers();
    virtual void refine(CrystalPtr target, RefinementType rType);


    static void setRotamerExponent(void *object, double exp)
    {
        if (exp < 0) exp = 0;
        Sidechain *side = static_cast<Sidechain *>(object);
        side->_exponent = exp;
        side->refreshRotamers();

    }

    static double getRotamerExponent(void *object)
    {
        return static_cast<Sidechain *>(object)->_exponent;
    }

protected:
    virtual bool shouldRefineMagicAxis(BondPtr bond);
    virtual AtomList topLevelAtoms()
    {
        return findAtoms("CB");
    }

    virtual std::string getClassName()
    {
        return "Sidechain";
    }

    virtual std::string getParserIdentifier()
    {
        return "side_" + i_to_str(_resNum);
    }

private:
    void refreshRotamers();
    bool _rotamerised;
    bool _canRefine;
    int _resNum;
    double _exponent; /* For rotamer weighting */
    PolymerWkr _myPolymer;
};

#endif /* defined(__vagabond__Sidechain__) */
