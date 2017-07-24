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

#define MAX_SCATTERING_DSTAR 3.0
#define ATOM_SAMPLING_COUNT (16)
#define PROTEIN_SAMPLING (1. / 3.)
#define WATER_RADIUS 0.6

class cFFTW3d;
typedef std::shared_ptr<cFFTW3d> FFTPtr;

class Crystal;
typedef std::shared_ptr<Crystal> CrystalPtr;

class Molecule;
typedef std::shared_ptr<Molecule> MoleculePtr;

class Monomer;
typedef std::shared_ptr<Monomer> MonomerPtr;
typedef std::weak_ptr<Monomer> MonomerWkr;

class Backbone;
typedef std::shared_ptr<Backbone> BackbonePtr;

class Sidechain;
typedef std::shared_ptr<Sidechain> SidechainPtr;

class Knotter;
typedef std::shared_ptr<Knotter> KnotterPtr;

class Polymer;
typedef std::shared_ptr<Polymer> PolymerPtr;
typedef std::weak_ptr<Polymer> PolymerWkr;

class Atom;
typedef std::shared_ptr<Atom> AtomPtr;
typedef std::weak_ptr<Atom> AtomWkr;

class Element;
typedef std::shared_ptr<Element> ElementPtr;

class Model;
class Absolute;
class Bond;

typedef std::shared_ptr<Absolute> AbsolutePtr;
typedef std::shared_ptr<Model> ModelPtr;
typedef std::shared_ptr<Bond> BondPtr;

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
