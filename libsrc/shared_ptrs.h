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
#define PROTEIN_SAMPLING (1. / 3.)
#define WATER_RADIUS 1.4

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

typedef enum
{
	MaskUnchecked = 0,
	MaskProtein = 1,
	MaskEmpty = 2,
} MaskType;

#endif
