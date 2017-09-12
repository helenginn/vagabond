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
#include <math.h>
#include <vector>

#define ANGLE_SAMPLING deg2rad(4.0)

#define MAX_SCATTERING_DSTAR 2.00
#define ATOM_SAMPLING_COUNT (36)
#define PROTEIN_SAMPLING (1. / 3.)
#define WATER_RADIUS 0.6

#define deg2rad(a) ((a) * M_PI / 180)
#define rad2deg(a) ((a) / M_PI * 180)

#define ToBondPtr(a) (std::static_pointer_cast<Bond>((a)))
#define ToAbsolutePtr(a) (std::static_pointer_cast<Absolute>((a)))
#define ToPolymerPtr(a) (std::static_pointer_cast<Polymer>((a)))

class FFT;
typedef std::shared_ptr<FFT> FFTPtr;

class Crystal;
typedef std::shared_ptr<Crystal> CrystalPtr;

class Molecule;
typedef std::shared_ptr<Molecule> MoleculePtr;

class Monomer;
typedef std::shared_ptr<Monomer> MonomerPtr;
typedef std::weak_ptr<Monomer> MonomerWkr;

class AtomGroup;
class Backbone;
typedef std::shared_ptr<Backbone> BackbonePtr;
typedef std::shared_ptr<AtomGroup> AtomGroupPtr;

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

class RefinementGridSearch;
class RefinementStepSearch;
class RefinementStrategy;
class RefinementSnake;
class NelderMead;
typedef std::shared_ptr<RefinementStepSearch> RefinementStepSearchPtr;
typedef std::shared_ptr<RefinementGridSearch> RefinementGridSearchPtr;
typedef std::shared_ptr<RefinementStrategy> RefinementStrategyPtr;
typedef std::shared_ptr<RefinementSnake> RefinementSnakePtr;
typedef std::shared_ptr<NelderMead> NelderMeadPtr;

class CSV;
class PNGFile;
class TextManager;
typedef std::shared_ptr<PNGFile> PNGFilePtr;
typedef std::shared_ptr<TextManager> TextManagerPtr;
typedef std::shared_ptr<CSV> CSVPtr;

class Sampler;
typedef std::shared_ptr<Sampler> SamplerPtr;

typedef std::vector<AtomWkr> AtomList;

typedef enum
{
	MaskUnchecked = 0,
	MaskProtein = 1,
	MaskEmpty = 2,

	MaskFree = 8,
	MaskWork = 9,
} MaskType;

#endif
