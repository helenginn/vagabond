//
//  AtomGroup.cpp
//  vagabond
//
//  Created by Helen Ginn on 25/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "AtomGroup.h"
#include "Atom.h"
#include "Element.h"
#include "Bond.h"
#include <sstream>
#include "Crystal.h"
#include <iomanip>
#include "maths.h"
#include "Shouter.h"
#include "../libccp4/ccp4_spg.h"
#include "FlexRegion.h"

AtomPtr AtomGroup::findAtom(std::string atomType)
{
    for (int i = 0; i < atomCount(); i++)
    {
        if (atom(i)->getAtomName() == atomType)
        {
            return atom(i);
        }
    }

    return AtomPtr();
}

AtomPtr AtomGroup::findAtom(std::string atomType, std::string confID)
{
    AtomList atoms = findAtoms(atomType);

    for (int i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].expired())
        {
            continue;
        }

        if (atoms[i].lock()->getAlternativeConformer() == confID)
        {
            return atoms[i].lock();
        }
    }

    return AtomPtr();
}

std::map<std::string, int> AtomGroup::conformerMap()
{
    std::map<std::string, int> conformerList;

    for (int i = 0; i < atomCount(); i++)
    {
        std::string conformer = atom(i)->getAlternativeConformer();

        if (!conformerList.count(conformer))
        {
            conformerList[conformer] = 0;
        }

        conformerList[conformer]++;
    }

    return conformerList;
}

int AtomGroup::conformerCount()
{
    std::map<std::string, int> conformerList = conformerMap();

    return conformerList.size();
}

std::string AtomGroup::conformer(int i)
{
    if (i > conformerMap().size()) return "";

    std::map<std::string, int> conformerList = conformerMap();
    std::map<std::string, int>::iterator it = conformerList.begin();

    for (int j = 0; j < i; j++) it++;

    return it->first;
}

AtomList AtomGroup::findAtoms(std::string atomType)
{
    AtomList list;

    for (int i = 0; i < atomCount(); i++)
    {
        if (atom(i)->getAtomName() == atomType)
        {
            list.push_back(atom(i));
        }
    }

    return list;
}

double AtomGroup::totalElectrons()
{
    double total = 0;

    for (int i = 0; i < atomCount(); i++)
    {
        total += atom(i)->getElement()->electronCount();
    }

    return total;
}

std::string AtomGroup::getPDBContribution(PDBType pdbType, CrystalPtr crystal)
{
    std::ostringstream stream;
    int numConf = 0;

    if (!atomCount())
    {
        return "";
    }

    if (pdbType == PDBTypeEnsemble)
    {
        /* Get the total number of conformers to worry about */
        std::vector<BondSample> *samples = atom(0)->getModel()->getManyPositions(BondSampleThorough);

        numConf = samples->size();

        for (int j = 0; j < numConf; j++)
        {
            stream << "MODEL " << std::setw(8) << j << std::setw(66) << " " << std::endl;

            for (int i = 0; i < atomCount(); i++)
            {
                if (!atom(i)->getMonomer())
                {
                    continue;
                }

                if (atom(i)->getWeighting() <= 0)
                {
                    continue;
                }

                stream << atom(i)->getPDBContribution(j);
            }
            
            stream << "TER" << std::setw(80) << " " << std::endl;
            stream << "ENDMDL" << std::setw(80) << " " << std::endl;
        }

        return stream.str();
    }

    for (int i = 0; i < atomCount(); i++)
    {
        bool samePos = (pdbType == PDBTypeSamePosition);
        bool sameB = (pdbType == PDBTypeSameBFactor);
        stream << atom(i)->averagePDBContribution(samePos, sameB);
        if (crystal)
        {
            stream << atom(i)->anisouPDBLine(crystal);
        }
    }

    return stream.str();
}

void AtomGroup::setUseAbsolute()
{
    for (int i = 0; i < atomCount(); i++)
    {
        atom(i)->setKeepModel();
    }
}

void AtomGroup::addAtomsFrom(AtomGroupPtr child)
{
    for (int i = 0; i < child->atomCount(); i++)
    {
        addAtom(child->atom(i));
    }
}

double AtomGroup::getAverageDisplacement()
{
    double sum = 0;
    double count = 0;

    for (int i = 0; i < atomCount(); i++)
    {
        if (atom(i)->getElement()->electronCount() <= 1)
        {
            continue;
        }

        double val = atom(i)->posDisplacement();

        sum += val;
        count++;
    }

    return sum / count;
}

double AtomGroup::getAverageBFactor(bool initial)
{
    double sum = 0;
    double count = 0;

    for (int i = 0; i < atomCount(); i++)
    {
        if (atom(i)->getElement()->electronCount() <= 1)
        {
            continue;
        }

        if (initial)
        {
            sum += atom(i)->getInitialBFactor();
            count++;
        }
        else
        {
            if (atom(i)->getModel()->isBond())
            {
                BondPtr bond = ToBondPtr(atom(i)->getModel());
                double val = bond->getMeanSquareDeviation();
                sum += val;
                count++;
            }
        }
    }

    return sum / count;
}

AtomGroup::AtomGroup()
{
    _beenTied = false;
    _timesRefined = 0;
}

void AtomGroup::propagateChange()
{
    for (int i = 0; i < atomCount(); i++)
    {
        atom(i)->getModel()->propagateChange(0);
    }
}

void AtomGroup::refreshPositions(bool quick)
{
    for (int i = 0; i < atomCount(); i++)
    {
        if (!atom(i)) continue;

        if (quick)
        {
            atom(i)->getModel()->propagateChange(0);
            atom(i)->getModel()->getFinalPositions();
            continue;
        }

        atom(i)->getModel()->propagateChange(-1, true);
    }
}


int AtomGroup::totalElectrons(int *fcWeighted)
{
    double sum = 0;
    double weighted = 0;

    for (int i = 0; i < atomCount(); i++)
    {
        double e = atom(i)->getElement()->electronCount();
        sum += e;
        double weight = atom(i)->getWeighting();
        weighted += e * weight;
    }

    *fcWeighted = weighted;

    return sum;
}

void AtomGroup::setWeighting(double value)
{
    for (int i = 0; i < atomCount(); i++)
    {
        atom(i)->setWeighting(value);
    }
}

void AtomGroup::resetMagicAxes()
{
}

AtomList AtomGroup::topLevelAtoms()
{
    if (!atomCount()) return AtomList();

    AtomPtr topAtom = atom(0);

    while (true)
    {
        if (!topAtom->getModel()->isBond())
        {
            break;
        }

        BondPtr bond = ToBondPtr(topAtom->getModel());

        if (!hasAtom(bond->getMajor()))
        {
            break;
        }

        topAtom = bond->getMajor();
    }

    AtomList list;
    list.push_back(topAtom);

    return list;
}

bool AtomGroup::hasAtom(AtomPtr anAtom)
{
    bool found = false;

    if (!anAtom) return false;

    for (int i = 0; i < atomCount(); i++)
    {
        if (atom(i) == anAtom)
        {
            found = true;
        }
    }

    return found;
}

void AtomGroup::setTargetRefinement(CrystalPtr target, RefinementType rType)
{
    _target = target;
    _rType = rType;
}

void AtomGroup::privateRefine()
{
    std::cout << "Refining one residue." << std::endl;
    refine(_target, _rType);
}

void AtomGroup::refine(CrystalPtr target, RefinementType rType)
{
    AtomList topAtoms = topLevelAtoms();
    bool refineAngles = shouldRefineAngles();
    _timesRefined++;

    ScoreType scoreType = ScoreTypeModelPos;
    int maxTries = 0;
    int bondNum = 4;
    double degrees = 0;

    switch (rType) {
        case RefinementModelPos:
        scoreType = ScoreTypeModelPos;
        maxTries = 60;
        degrees = 4;
        break;

        case RefinementFine:
        scoreType = ScoreTypeCorrel;
        maxTries = 6;
        degrees = 4;
        bondNum = 2;
        refineAngles = false;
        break;

        case RefinementFlexibility:
        scoreType = ScoreTypeModelPos;
        maxTries = 3;
        degrees = 4;
        break;

        default:
        shout_at_helen("Unimplemented refinement option?");
        break;
    }

    if (refineAngles)
    {
        bondNum = 3;
    }

    for (int n = 0; n < topAtoms.size(); n++)
    {
        AtomPtr topAtom = topAtoms[n].lock();

        if (n == 1) std::cout << "'" << std::flush;

        while (hasAtom(topAtom))
        {
            if (!topAtom->getModel()->isBond())
            {
                break;
            }

            BondPtr bond = ToBondPtr(topAtom->getModel());

            int groups = bond->downstreamAtomGroupCount();

            if (!groups)
            {
                break;
            }

            if (!bond->isRefinable())
            {
                break;
            }

            for (int k = 0; k < 1; k++)
            {
                if (shouldRefineMagicAxis(bond))
                {
                    bond->calculateMagicAxis();
                }

                int count = 0;

                BondPtr topBond;

                while (rType == RefinementFlexibility && count < maxTries)
                {
                    FlexRegion flexer;
                    flexer.setup();
                    flexer.addBond(ToBondPtr(bond), 6);
                    flexer.addSingleBondParameters();
                    flexer.sample();
                    count++;
                }

                if (rType == RefinementFlexibility)
                {
                    rType = RefinementModelPos;
                    count = 0;
                    maxTries = 60;
                }

                bool changed = true;
                bool addFlex = (rType == RefinementFine);

                while (changed && count < maxTries)
                {
                    bond->setActiveGroup(k);
                    setupNelderMead();
                    setCrystal(target);
                    setCycles(16);
                    topBond = setupTorsionSet(bond, k, bondNum,
                                              deg2rad(degrees), deg2rad(0.04),
                                              refineAngles, addFlex);
                    setScoreType(scoreType);

                    if (rType == RefinementModelPos)
                    {
                        setSilent();
                    }

                    setJobName("torsion_" +  bond->shortDesc());
                    changed = sample();
                    count++;
                }

                if (!topBond)
                {
                    topAtom = AtomPtr();
                    continue;
                }

                topAtom = topBond->getMinor();
            }
        }
    }

    refreshPositions(true);
}

double AtomGroup::scoreWithMap(ScoreType scoreType, CrystalPtr crystal)
{
    FFTPtr map = crystal->getFFT();
    mat3x3 real2Frac = crystal->getReal2Frac();

    return scoreWithMap(scoreType, map, real2Frac);
}

double AtomGroup::scoreWithMap(ScoreType scoreType, FFTPtr map, mat3x3 real2Frac)
{
    return scoreWithMap(_atoms, scoreType, map, real2Frac);
}

double AtomGroup::scoreWithMap(std::vector<AtomPtr> atoms, ScoreType scoreType,
                               FFTPtr map, mat3x3 real2Frac)
{
    atoms[0]->getModel()->getDistribution(true);
    vec3 zero = atoms[0]->getModel()->getAbsolutePosition();
    double maxDistance = 0;

    for (int i = 1; i < atoms.size(); i++)
    {
        /* Refresh absolute position */
        atoms[i]->getModel()->getDistribution();
        vec3 offset = atoms[i]->getModel()->getAbsolutePosition();

        vec3 diff = vec3_subtract_vec3(offset, zero);
        double distance = vec3_length(diff);

        if (distance > maxDistance)
        {
            maxDistance = distance;
        }
    }

    double scales = 0.6;
    double n = 2 * (maxDistance + 3.0) / scales;
//    n = 24;

    FFTPtr segment = FFTPtr(new FFT());
    segment->create(n + 0.5);
    segment->setScales(scales);
    mat3x3 basis = make_mat3x3();
    double toReal = 1 / (scales*n);
    mat3x3_scale(&basis, toReal, toReal, toReal);
    segment->createFFTWplan(1);

    for (int i = 0; i < atoms.size(); i++)
    {
        atoms[i]->addToMap(segment, basis, zero);
    }

//    segment->normalise();
    mat3x3_mult_vec(real2Frac, &zero);

    std::vector<double> xs, ys;

    double cutoff = FFT::score(map, segment, zero, &xs, &ys);

    FFTPtr obsSeg = FFTPtr(new FFT(*segment));
    obsSeg->setAll(0);

    /* n.b. this is fucked. please unfuck before continuing. */
//    FFT::score(map, obsSeg, zero, NULL, NULL, MapScoreTypeCopyToSmaller);

    /*

    segment->shiftToCentre();
    segment->normalise();
    segment->printSlice();

    segment->fft(-1);
    segment->writeReciprocalToFile("segment_calc.mtz");

    obsSeg->shiftToCentre();
    obsSeg->normalise();
    obsSeg->fft(-1);
    obsSeg->writeReciprocalToFile("segment_obs.mtz");
*/

    /*
    double correl = correlation(xs, ys, cutoff);
    std::cout << "Correlation: " << correl << std::endl;
    std::cout << "Cutoff: " << cutoff << std::endl;

    if (correl > 0.65)
    {
        for (int i = 0; i < xs.size(); i++)
        {
            if (ys[i] > cutoff)
            {
                std::cout << xs[i] << ", " << ys[i] << std::endl;
            }
        }

        exit(0);
    }
     */

    if (scoreType == ScoreTypeCorrel)
    {
        double correl = correlation(xs, ys, cutoff);
        return -correl;
    }
    else if (scoreType == ScoreTypeRFactor)
    {
        double rFactor = scaled_r_factor(xs, ys, cutoff);
        return rFactor;
    }
    else if (scoreType == ScoreTypeMultiply)
    {
        double mult = weightedMapScore(xs, ys);
        return -mult;
    }

    return 0;
}

void AtomGroup::addProperties()
{

}
