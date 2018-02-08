//
//  AtomGroup.h
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__AtomGroup__
#define __vagabond__AtomGroup__

#include <stdio.h>
#include <string>
#include <vector>
#include "shared_ptrs.h"
#include "Sampler.h"
#include "Parser.h"
#include <map>

class AtomGroup : public boost::enable_shared_from_this<AtomGroup>, public Sampler, public Parser
{
public:
    AtomPtr findAtom(std::string atomType);
    AtomPtr findAtom(std::string atomType, std::string confID);
    AtomList findAtoms(std::string atomType);

    double scoreWithMap(ScoreType scoreType, CrystalPtr crystal);
    static double scoreWithMap(std::vector<AtomPtr> atoms, ScoreType scoreType,
                               FFTPtr map, mat3x3 real2Frac);
    double scoreWithMap(ScoreType scoreType, FFTPtr map, mat3x3 real2Frac);

    void setMonomer(MonomerPtr monomer)
    {
        _monomer = monomer;
    }

    MonomerPtr getMonomer()
    {
        return _monomer.lock();
    }

    virtual void addAtom(AtomPtr atom)
    {
        std::vector<AtomPtr>::iterator it;
        it = std::find(_atoms.begin(), _atoms.end(), atom);

        if (it == _atoms.end())
        {
            _atoms.push_back(atom);
        }
    }

    long atomCount()
    {
        return _atoms.size();
    }

    bool hasAtom(AtomPtr anAtom);
    
        AtomPtr atom(int i)
    {
        return _atoms[i];
    }

    void addBond(BondPtr bond)
    {
        _bonds.push_back(bond);
    }

    int bondCount()
    {
        return _bonds.size();
    }

    BondPtr bond(int i)
    {
        return _bonds[i].lock();
    }

    double totalElectrons();
    double getAverageBFactor(bool initial = false);
    double getAverageDisplacement();

    std::string getPDBContribution(PDBType pdbType,
                                   CrystalPtr crystal = CrystalPtr());

    void setTied()
    {
        _beenTied = true;
    }

    virtual bool shouldRefineAngles()
    {
        return false;
    }

    void setUseAbsolute();
    int totalElectrons(int *fcWeighted);


    static double refine(void *object)
    {
        static_cast<AtomGroup *>(object)->privateRefine();
    }

    void setTargetRefinement(CrystalPtr target, RefinementType rType);
    virtual void refine(CrystalPtr target, RefinementType rType);
    void setWeighting(double value);
    void resetMagicAxes();
    int conformerCount();
    std::string conformer(int i);
    void propagateChange();
    void refreshPositions(bool quick = true);
protected:
    AtomGroup();
    void addAtomsFrom(AtomGroupPtr child);
    virtual AtomList topLevelAtoms();
    int _timesRefined;

    bool isTied()
    {
        return _beenTied;
    }

    virtual std::string getClassName()
    {
        return "AtomGroup";
    }

    virtual std::string getParserIdentifier()
    {
        return "AtomGroupSomething";
    }

    virtual void addProperties();

private:
    MonomerWkr _monomer;

    std::vector<BondWkr> _bonds;
    std::vector<AtomPtr> _atoms;

    bool _beenTied;
    CrystalPtr _target;
    RefinementType _rType;
   
    void privateRefine(); 
    std::map<std::string, int> conformerMap();

};

#endif /* defined(__vagabond__AtomGroup__) */
