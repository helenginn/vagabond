//
//  Bond.h
//  vagabond
//
//  Created by Helen Ginn on 23/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Bond__
#define __vagabond__Bond__

#include <stdio.h>
#include <string>
#include <vector>
#include "mat3x3.h"
#include "Distributor.h"
#include <iostream>
#include "Atom.h"
#include "Sampler.h"
#include "Model.h"
#include <mutex>

typedef struct
{
    AtomWkr atom;
    double geomRatio;
    double expectedAngle;
    double circlePortion;
    std::string *placeholder;
} AtomValue;

typedef struct
{
    std::vector<AtomValue> atoms;
    double torsionAngle;
    double torsionBlur;
    double magicPhi;
    double magicPsi;
    std::vector<BondSample> storedSamples;
    std::vector<AtomWkr> extraTorsionSamples;
} BondGroup;

#define INITIAL_KICK 0.01
#define INITIAL_DAMPENING 0.08

class Anisotropicator;

class Bond : public Model, public Sampler
{
public:
    Bond(AtomPtr major, AtomPtr minor, int group = 0);
    Bond(Bond &other);
    Bond();
    void initialize();
    void activate(AtomGroupPtr group = AtomGroupPtr(),
                  AtomPtr inherit = AtomPtr());
    void setupSampling();
    std::vector<AtomPtr> importantAtoms();
    ModelPtr reverse(BondPtr upstreamBond);
    void reverseDownstreamAtoms(int group);
    void resetBondDirection();
    bool isRefinable();
    void calculateMagicAxis();
    double magicAxisScore();
    static double magicAxisStaticScore(void *object);
    bool test();
    double getEffectiveOccupancy();

    AtomPtr getMajor()
    {
        return _major.lock();
    }

    AtomPtr getMinor()
    {
        return _minor.lock();
    }

    void setMajor(AtomPtr newMajor)
    {
        _major = newMajor;
    }

    void setMinor(AtomPtr newMinor);

    AtomPtr getHeavyAlign()
    {
        return _heavyAlign.lock();
    }

    AtomPtr getLightAlign()
    {
        return _lightAlign.lock();
    }

    static double getBondLength(void *object)
    {
        return static_cast<Bond *>(object)->_bondLength;
    }

    static void setBondLength(void *object, double length)
    {
        static_cast<Bond *>(object)->_bondLength = length;
    }

    double getTorsion(int group)
    {
        return _bondGroups[group].torsionAngle;
    }

    void setTorsionAtoms(AtomPtr heavyAlign = AtomPtr(),
                         AtomPtr lightAlign = AtomPtr(),
                         int groupNum = 0);
    virtual FFTPtr getDistribution(bool quick = false);

    virtual std::string getClassName()
    {
        return "Bond";
    }

    static void setTorsionBlur(void *object, double value)
    {
        Bond *bond = static_cast<Bond *>(object);
        bond->_bondGroups[bond->_activeGroup].torsionBlur = value;
        static_cast<Bond *>(object)->propagateChange(20);
    }

    static double getTorsionBlur(void *object)
    {
        Bond *bond = static_cast<Bond *>(object);
        return bond->_bondGroups[bond->_activeGroup].torsionBlur;
    }

    static double getTorsion(void *object)
    {
        Bond *bond = static_cast<Bond *>(object);
        return bond->_bondGroups[bond->_activeGroup].torsionAngle;
    }

    static void setTorsion(void *object, double value)
    {
        Bond *bond = static_cast<Bond *>(object);
        bond->_bondGroups[bond->_activeGroup].torsionAngle = value;
        static_cast<Bond *>(object)->propagateChange(20);
    }

    static double getMagicPsi(void *object)
    {
        Bond *bond = static_cast<Bond *>(object);
        return bond->_bondGroups[bond->_activeGroup].magicPsi;
    }

    static void setMagicPsi(void *object, double angle)
    {
        Bond *bond = static_cast<Bond *>(object);
        bond->_bondGroups[bond->_activeGroup].magicPsi = angle;
    }

    static double getMagicPhi(void *object)
    {
        Bond *bond = static_cast<Bond *>(object);
        return bond->_bondGroups[bond->_activeGroup].magicPhi;
    }

    static void setMagicPhi(void *object, double angle)
    {
        Bond *bond = static_cast<Bond *>(object);
        bond->_bondGroups[bond->_activeGroup].magicPhi = angle;
        }

    static double getOccupancy(void *object)
    {
        Bond *bond = static_cast<Bond *>(object);
        return bond->_occupancy;
    }

    static void setOccupancy(void *object, double value);

    static double getDampening(void *object)
    {
        return static_cast<Bond *>(object)->_dampening;
    }

    static void setDampening(void *object, double value)
    {
        Bond *bond = static_cast<Bond *>(object);
        bond->_dampening = value;
        bond->propagateChange(20);
    }

    static double getBendAngle(void *object);
    static void setBendAngle(void *object, double value);

    void addDownstreamAtom(AtomPtr atom, int group, bool skipGeometry = false);

    int downstreamAtomCount(int group)
    {
        return _bondGroups[group].atoms.size();
    }

    int downstreamAtomGroupCount()
    {
        return _bondGroups.size();
    }

    AtomPtr downstreamAtom(int group, int i)
    {
        return _bondGroups[group].atoms[i].atom.lock();
    }

    int downstreamAtomNum(AtomPtr atom, int *group)
    {
        for (int j = 0; j < downstreamAtomGroupCount(); j++)
        {
            for (int i = 0; i < downstreamAtomCount(j); i++)
            {
                if (downstreamAtom(j, i) == atom)
                {
                    *group = j;
                    return i;
                }
            }
        }

        *group = -1;
        return -1;
    }

    void setUsingTorsion(bool use)
    {
        _usingTorsion = use;
    }

    bool isUsingTorsion()
    {
        return _usingTorsion;
    }

    bool isNotJustForHydrogens();

    double getExpectedAngle();
    double getGeomRatio(int n, int i)
    {
        return _bondGroups[n].atoms[i].geomRatio;
    }

    double getCirclePortion(int n, int i)
    {
        return _bondGroups[n].atoms[i].circlePortion;
    }

    static double getCirclePortion(void *object);
    static void setCirclePortion(void *object, double value);

    void setGeomRatio(int n, int i, double value)
    {
        _bondGroups[n].atoms[i].geomRatio = value;
    }

    void setActiveGroup(int newGroup)
    {
        if (_activeGroup == newGroup)
        {
            return;
        }

        _activeGroup = newGroup;
        propagateChange(0);
    }

    int getActiveGroup()
    {
        return _activeGroup;
    }

    BondGroup *getBondGroup(int i)
    {
        return &_bondGroups[i];
    }

    void setOccupancyMult(double mult)
    {
        _occMult = mult;
        propagateChange();
    }

    std::string description();
    std::string shortDesc();
    std::string getPDBContribution();
    ModelPtr getParentModel();
    bool splitBond(int start = 0);

    void setFixed(bool fixed)
    {
        _fixed = fixed;
    }

    bool isFixed()
    {
        return _fixed;
    }

    void addExtraTorsionSample(AtomPtr atom, int group)
    {
        _bondGroups[group].extraTorsionSamples.push_back(atom);
    }

    int extraTorsionSampleCount(int group)
    {
        return _bondGroups[group].extraTorsionSamples.size();
    }

    virtual double getMeanSquareDeviation();

    double getFlexibilityPotential();

    AtomPtr extraTorsionSample(int group, int i)
    {
        return _bondGroups[group].extraTorsionSamples[i].lock();
    }

    vec3 getAbsolutePosition()
    {
        return _absolute;
    }

    void setAnchored()
    {
        _anchored = true;
        _changedPos = false;
        _changedSamples = false;
    }


    /* Will aim to define torsion basis as:
     x: along line of 0º torsion angle.
     y: completes the right-handed coordinate system
     z: along bond direction, from heavy-to-light alignment atoms.
     */
    mat3x3 makeTorsionBasis(vec3 hPos, vec3 maPos,
                            vec3 miPos, vec3 lPos, double *newAngle = NULL);

    void recalculateTorsion(AtomPtr heavy, double value);

    virtual void propagateChange(int depth = -1, bool refresh = false);
    std::vector<BondSample> *getManyPositions();

    std::vector<vec3> polymerCorrectedPositions();

    static void useMutex()
    {
        _useMutex = true;
    }

    double getMultOccupancy()
    {
        return _occupancy * _occMult;
    }

    void setRefineBondAngle(bool value = true)
    {
        _refineBondAngle = value;
    }

    bool getRefineBondAngle()
    {
        return _refineBondAngle * !isFixed();
    }

    static void encodeBondGroup(void *bond, void *bondGroup,
                                std::ofstream &stream, int indent);
    static char *decodeBondGroup(void *bond, void *bondGroup, char *block);
protected:

    AtomWkr _minor;

    virtual std::string getParserIdentifier()
    {
        return "bond_" + shortDesc();
    }

    virtual void addProperties();
    virtual void linkReference(ParserPtr object, std::string category);
    virtual void postParseTidy();    
private:
    std::string _shortDesc;

    AtomWkr _major;

    AtomWkr _heavyAlign;
    AtomWkr _lightAlign;

    double _bondLength;

    /* Downstream groups of bonds */
    std::vector<BondGroup> _bondGroups;
    std::vector<AtomWkr> _magicAxisAtoms;

    int magicAtomCount()
    {
        return _magicAxisAtoms.size();
    }

    AtomPtr getMagicAtom(int i)
    {
        return _magicAxisAtoms[i].lock();
    }

    double _dampening;
    bool _activated;
    int _activeGroup;
    double _blurTotal;
    double _occupancy;
    double _occMult;

    /* Should not be refined */
    bool _fixed;

    /* Had a non-NULL atom input as major or minor */
    bool _disabled;

    /* Has been set as an anchor, will not respond to 'propagate change'*/
    bool _anchored;

    /* Grab bond length from the atom types of major/minor */
    void deriveBondLength();
    /* And a given bond angle for downstream atom n */
    double deriveBondAngle(AtomPtr atom);

    /* Bond direction only used when a torsion angle can't be
     * calculated because it's connected to an Absolute PDB.
     * Otherwise use as a reference for torsion matrix updates. */
    vec3 _bondDirection;

    vec3 positionFromTorsion(mat3x3 torsionBasis, double angle,
                             double ratio, vec3 start);
    std::vector<BondSample> getCorrectedAngles(std::vector<BondSample> *prevs,
                                               double circleAdd,
                                               double myTorsion, double ratio);

    BondPtr duplicateDownstream(BondPtr newBranch, int groupNum, int start);
    bool _usingTorsion;

    /* Flag to say whether recalculation should occur */
    bool _changedPos, _changedSamples;
    bool _refineBondAngle;

    mat3x3 getMagicMat(vec3 direction);


};

#endif /* defined(__vagabond__Bond__) */
