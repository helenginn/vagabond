//
//  Molecule.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Molecule__
#define __vagabond__Molecule__

#include <stdio.h>
#include <vector>
#include "mat3x3.h"
#include "shared_ptrs.h"
#include "Sampler.h"
#include "AtomGroup.h"
#include <string>
#include <iostream>

struct vec3;

class Molecule : public AtomGroup
{
public:
    Molecule();

    void addToMap(FFTPtr fft, mat3x3 _real2frac);

    virtual void summary();
    virtual void tieAtomsUp();
    virtual void refine(CrystalPtr target, RefinementType rType);
    void tiedUpScattering(double *tied, double *all);
    virtual std::string makePDB(PDBType pdbType, CrystalPtr crystal);
    virtual void graph(std::string graphName) {};
    virtual void differenceGraphs(std::string graphName, CrystalPtr diffCryst) {};
    void resetInitialPositions();
    void setAnchors();
    void makePowderList();

    virtual void reportParameters();

    double getAbsoluteBFacSubtract()
    {
            return _absoluteBFacSubtract;
    }

    void setAbsoluteBFacSubtract(double subtract)
    {
         _absoluteBFacSubtract = subtract;
    }

    double getAbsoluteBFacMult()
    {
          return _absoluteBFacMult;
    }

    void setAbsoluteBFacMult(double mult)
    {
         _absoluteBFacMult = mult;
    }

    void setChainID(std::string chain)
    {
        _chainID = chain;
    }

    std::string getChainID()
    {
        return _chainID;
    }

    std::vector<vec3> getCentroidOffsets()
    {
        return _centroidOffsets;
    }

    std::vector<mat3x3> getRotationCorrections()
    {
        return _rotations;
    }

    std::vector<vec3> getRotationCentres()
    {
        return _centroids;
    }

    std::vector<vec3> getTransTensorOffsets()
    {
        return _transTensorOffsets;
    }

    std::vector<mat3x3> getExtraRotations();

    virtual std::string getClassName()
    {
        return "Molecule";
    }

    bool isPolymer()
    {
        return (getClassName() == "Polymer");
    }

    MoleculePtr shared_from_this()
    {
        AtomGroupPtr groupPtr = AtomGroup::shared_from_this();
        return ToMoleculePtr(groupPtr);
    }
protected:
    virtual std::string getParserIdentifier()
    {
        return "chain_" + _chainID; 
    }

    virtual void addProperties();
    virtual void postParseTidy();    

    std::vector<vec3> _centroidOffsets;
    std::vector<vec3> _centroids; // after offset correction
    std::vector<mat3x3> _rotations;

    std::vector<vec3> _transTensorOffsets;
    std::vector<mat3x3> _extraRotationMats; // currently unused

    virtual void calculateExtraRotations() {};

    // this axis calculates the angular response to the reaction sphere
    vec3 _magicRotAxis;

    // this axis is that of the rotation matrices applied to the structure
    vec3 _rotationAxis;    

    double _rotationAngle;

    void setChangedRotation()
    {
        _changedRotations = true;
    }
    
private:
    double _absoluteBFacSubtract;
    double _absoluteBFacMult;

    bool _changedRotations;
    std::string _chainID;
};

#endif /* defined(__vagabond__Molecule__) */
