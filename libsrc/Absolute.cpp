//
//  Absolute.cpp
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Absolute.h"
#include "shared_ptrs.h"
#include "Atom.h"
#include <math.h>
#include "fftw3d.h"
#include "Element.h"
#include <iostream>
#include "Monomer.h"
#include "Polymer.h"
#include "maths.h"
#include "Crystal.h"
#include "Anisotropicator.h"
#include <iomanip>
#include <iostream>

Absolute::Absolute()
{
    _position = make_vec3(0, 0, 0);
    _bFactor = 0;
    _element = "";
    _occupancy = 1;
    _hetatm = false;
    _usingTensor = false;
    _tensor = make_mat3x3();
    _isOfManyPositions = false;
}

mat3x3 Absolute::getRealSpaceTensor()
{
    if (_isOfManyPositions)
    {
        getAnisotropy(true);
        return Model::getRealSpaceTensor();
    }

    return _tensor;
}

void Absolute::getAnisotropy(bool withKabsch)
{
    if (_isOfManyPositions)
    {
        Model::getAnisotropy(true);
    }
}

Absolute::Absolute(vec3 pos, double bFac, std::string element, double occValue)
{
    _position = pos;
    _bFactor = bFac;
    _element = element;
    trim(_element);
    _occupancy = occValue;
    _hetatm = false;
    _usingTensor = false;
    _tensor = make_mat3x3();
    _isOfManyPositions = false;
}

void Absolute::makeAtom()
{
    AtomPtr myAtom = AtomPtr(new Atom());
    myAtom->setModel(shared_from_this());
    myAtom->setInitialPosition(_position);
    myAtom->setInitialBFactor(_bFactor);
    myAtom->setPDBPosition(_position);
    myAtom->setAtomNum(_atomNum);
    myAtom->setAlternativeConformer(_conformer);
    ElementPtr element = Element::getElement(_element);
    myAtom->setElement(element);
    myAtom->setAtomName(_atomName);
    myAtom->findAtomType(_resName);
    myAtom->setOriginalOccupancy(_occupancy);

    _atom = myAtom;
}

void Absolute::addToMolecule(MoleculePtr molecule)
{
    makeAtom();
    molecule->addAtom(_atom);
    setMolecule(molecule);
}

double Absolute::getExpValue(void *object, double x, double y, double z)
{
    Absolute *me = static_cast<Absolute *>(object);
    double aniso = 0;
    double mult = 1;
    
    if (me->hasMolecule())
    {
        MoleculePtr molecule = me->getMolecule();
        mult = molecule->getAbsoluteBFacMult();
    }

    if (me->_usingTensor)
    {
        mat3x3 scaledTensor = me->getRealSpaceTensor();
        mat3x3_mult_scalar(&scaledTensor, mult);
        vec3 recipVec = make_vec3(x, y, z);
        mat3x3_mult_vec(scaledTensor, &recipVec);
        recipVec.x *= x;
        recipVec.y *= y;
        recipVec.z *= z;
        double multByTranspose = recipVec.x + recipVec.y + recipVec.z;
        aniso = exp((2 * M_PI * M_PI) * -(multByTranspose));

        return aniso * me->_occupancy;
    }

    double distSq =    (x * x + y * y + z * z);

    double bf = me->_bFactor;

    if (me->hasMolecule())
    {
        MoleculePtr molecule = me->getMolecule();
        double subtract = molecule->getAbsoluteBFacSubtract();
        bf -= subtract;
        bf *= mult;
    }

    double exponent = (-0.25) * bf * distSq;
    double value = exp(exponent);
    value *= me->_occupancy;

    return value;
}

/*  Absolute distribution only needs to be the blurring due to the atomic
 *  B factor. The position should be provided by a different function. */

FFTPtr Absolute::getDistribution(bool quick)
{
    double n = ATOM_SAMPLING_COUNT;
    double scale = 2 * MAX_SCATTERING_DSTAR;

    prepareDistribution(n, scale, this, Absolute::getExpValue);

    return getDistributionCopy();
}

std::vector<BondSample> *Absolute::getManyPositions(BondSampleStyle style)
{
    std::vector<BondSample> *bondSamples = &_bondSamples;
    bondSamples->clear();

    if (style == BondSampleStatic)
    {
        BondSample sample;
        sample.basis = make_mat3x3();
        sample.occupancy = 1;
        sample.torsion = 0;
        sample.old_start = make_vec3(0, 0, 0);
        sample.start = _position;
        bondSamples->push_back(sample);
        return bondSamples;
    }

        /* B factor isotropic only atm, get mean square displacement in
     * each dimension. */
    double meanSqDisp = getBFactor() / (8 * M_PI * M_PI);
    meanSqDisp = sqrt(meanSqDisp);

    double occTotal = 0;

    int samples = 101;
    int rnd = 1;

    std::vector<vec3> points;
    double offset = 2. / (double)samples;
    double increment = M_PI * (3.0 - sqrt(5));

    _sphereAngles.clear();

    double m = meanSqDisp;

    for (int i = 0; i < samples; i++)
    {
        double y = (((double)i * offset) - 1) + (offset / 2);
        double r = sqrt(1 - y * y);

        double phi = (double)((i + rnd) % samples) * increment;

        double x = cos(phi) * r;
        double z = sin(phi) * r;

        vec3 point = make_vec3(x * m, y * m, z * m);

        points.push_back(point);
        _sphereAngles.push_back(point);
    }

    for (int i = 0; i < points.size(); i++)
    {
        vec3 full = vec3_add_vec3(points[i], _position);
        double occ = 1;
        occTotal += occ;

        BondSample sample;
        sample.basis = make_mat3x3();
        sample.occupancy = occ;
        sample.torsion = 0;
        sample.old_start = make_vec3(0, 0, 0);
        sample.start = full;
        bondSamples->push_back(sample);
    }

    for (int i = 0; i < bondSamples->size(); i++)
    {
        (*bondSamples)[i].occupancy /= occTotal;
    }

    return bondSamples;
}

void Absolute::addToMonomer(MonomerPtr monomer)
{
    makeAtom();

    monomer->addAtom(_atom);
    monomer->getPolymer()->addAtom(_atom);
    setMolecule(monomer->getPolymer());

    Model::addToMonomer(monomer);
}

double Absolute::getMeanSquareDeviation()
{
    return _bFactor;
}

vec3 Absolute::getStaticPosition()
{
    return _position;
}

void Absolute::setTensor(mat3x3 tensor, CrystalPtr crystal)
{
    _tensor = tensor;
    _usingTensor = true;

//    mat3x3 toNormal = crystal->getHKL2Real();
//    _realSpaceTensor = mat3x3_mult_mat3x3(toNormal, _tensor);
    _realSpaceTensor = _tensor;

    Anisotropicator tropicator;
    tropicator.setTensor(_realSpaceTensor);
    vec3 longestAxis = tropicator.longestAxis();

    getAtom()->setEllipsoidLongestAxis(longestAxis);
    getAtom()->setTensor(_tensor);
}

void Absolute::addProperties()
{
    addDoubleProperty("occupancy", &_occupancy);
    addStringProperty("chain", &_chainID);
    addStringProperty("res_name", &_resName);
    addStringProperty("atom_name", &_atomName);
    addStringProperty("conformer", &_conformer);
    addIntProperty("res_num", &_resNum);
    addIntProperty("atom_num", &_atomNum);
    addBoolProperty("hetatm", &_hetatm);
    addBoolProperty("using_tensor", &_usingTensor);
    addVec3Property("position", &_position);
    addDoubleProperty("bfactor", &_bFactor);
    addBoolProperty("many_positions", &_isOfManyPositions);

    // missing tensor and sampled absolute positions from manyPositions
}

