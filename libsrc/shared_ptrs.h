//
//  shared_ptrs.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef vagabond_shared_ptrs_h
#define vagabond_shared_ptrs_h

#include <memory>

#define ATOM_MAX_RADIUS (2.0)
#define ATOM_SAMPLING (1. / 4.)

#define MAX_SCATTERING_DSTAR 6
#define ATOM_SAMPLING_DSTAR (1. / 4.)
#define ATOM_SAMPLING_COUNT (16)
#define PROTEIN_SAMPLING (1. / 4.)
#define WATER_RADIUS 0.6

class cFFTW3d;
typedef std::shared_ptr<cFFTW3d> FFTPtr;

class Crystal;
typedef std::shared_ptr<Crystal> CrystalPtr;

class Molecule;
typedef std::shared_ptr<Molecule> MoleculePtr;

class Atom;
typedef std::shared_ptr<Atom> AtomPtr;

class Element;
typedef std::shared_ptr<Element> ElementPtr;

class Model;
class Absolute;

typedef std::shared_ptr<Absolute> AbsolutePtr;
typedef std::shared_ptr<Model> ModelPtr;

class Bucket;
class BucketUniform;

typedef std::shared_ptr<Bucket> BucketPtr;
typedef std::shared_ptr<BucketUniform> BucketUniformPtr;

class Dataset;
class Diffraction;
class DiffractionMtz;
typedef std::shared_ptr<Dataset> DatasetPtr;
typedef std::shared_ptr<Diffraction> DiffractionPtr;
typedef std::shared_ptr<DiffractionMtz> DiffractionMtzPtr;

class Object;
typedef std::shared_ptr<Object> ObjectPtr;

typedef enum
{
	MaskUnchecked = 0,
	MaskProtein = 1,
	MaskEmpty = 2,
} MaskType;

#endif
