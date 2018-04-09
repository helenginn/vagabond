//
//  shared_ptrs.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef vagabond_shared_ptrs_h
#define vagabond_shared_ptrs_h

#define BOOST_DISABLE_ASSERTS
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <math.h>
#include <vector>

#define ANGLE_SAMPLING deg2rad(4.0)

#define MAX_SCATTERING_DSTAR 1.1
#define ATOM_SAMPLING_COUNT (24)
#define PROTEIN_SAMPLING (0.5)
#define WATER_RADIUS 0.6
#define FUTURE_MAGIC_ATOMS 12

#define deg2rad(a) ((a) * M_PI / 180)
#define rad2deg(a) ((a) / M_PI * 180)

#define ToBondPtr(a) (boost::static_pointer_cast<Bond>((a)))
#define ToAbsolutePtr(a) (boost::static_pointer_cast<Absolute>((a)))
#define ToModelPtr(a) (boost::static_pointer_cast<Model>((a)))
#define ToPolymerPtr(a) (boost::static_pointer_cast<Polymer>((a)))
#define ToMoleculePtr(a) (boost::static_pointer_cast<Molecule>((a)))
#define ToSidechainPtr(a) (boost::static_pointer_cast<Sidechain>((a)))
#define ToBackbonePtr(a) (boost::static_pointer_cast<Backbone>((a)))
#define ToCrystalPtr(a) (boost::static_pointer_cast<Crystal>((a)));
#define ToAtomPtr(a) (boost::static_pointer_cast<Atom>((a)));
#define ToAtomGroupPtr(a) (boost::static_pointer_cast<AtomGroup>((a)));
#define ToMonomerPtr(a) (boost::static_pointer_cast<Monomer>((a)));
#define ToParserPtr(a) (boost::static_pointer_cast<Parser>((a)));

#define ToThingPtr(a) (boost::static_pointer_cast<Thing>((a)));

class Notifiable;
typedef boost::shared_ptr<Notifiable> NotifiablePtr;

class Options;
typedef boost::shared_ptr<Options> OptionsPtr;

class FFT;
typedef boost::shared_ptr<FFT> FFTPtr;

class Crystal;
typedef boost::shared_ptr<Crystal> CrystalPtr;
typedef boost::weak_ptr<Crystal> CrystalWkr;

class Molecule;
typedef boost::shared_ptr<Molecule> MoleculePtr;
typedef boost::weak_ptr<Molecule> MoleculeWkr;

class Monomer;
typedef boost::shared_ptr<Monomer> MonomerPtr;
typedef boost::weak_ptr<Monomer> MonomerWkr;

class AtomGroup;
class Backbone;
typedef boost::shared_ptr<Backbone> BackbonePtr;
typedef boost::shared_ptr<AtomGroup> AtomGroupPtr;

class Sidechain;
typedef boost::shared_ptr<Sidechain> SidechainPtr;

class Knotter;
typedef boost::shared_ptr<Knotter> KnotterPtr;

class Polymer;
typedef boost::shared_ptr<Polymer> PolymerPtr;
typedef boost::weak_ptr<Polymer> PolymerWkr;

class Atom;
typedef boost::shared_ptr<Atom> AtomPtr;
typedef boost::weak_ptr<Atom> AtomWkr;

class Element;
typedef boost::shared_ptr<Element> ElementPtr;

class Model;
class Absolute;
class Bond;

typedef boost::shared_ptr<Absolute> AbsolutePtr;
typedef boost::shared_ptr<Model> ModelPtr;
typedef boost::shared_ptr<Bond> BondPtr;
typedef boost::weak_ptr<Bond> BondWkr;

class Bucket;
class BucketBulkSolvent;

typedef boost::shared_ptr<Bucket> BucketPtr;
typedef boost::shared_ptr<BucketBulkSolvent> BucketBulkSolventPtr;

class Dataset;
class Diffraction;
class DiffractionMtz;
typedef boost::shared_ptr<Dataset> DatasetPtr;
typedef boost::shared_ptr<Diffraction> DiffractionPtr;
typedef boost::shared_ptr<DiffractionMtz> DiffractionMtzPtr;

class Object;
typedef boost::shared_ptr<Object> ObjectPtr;

class RefinementGridSearch;
class RefinementStepSearch;
class RefinementStrategy;
class RefinementSnake;
class NelderMead;
typedef boost::shared_ptr<RefinementStepSearch> RefinementStepSearchPtr;
typedef boost::shared_ptr<RefinementGridSearch> RefinementGridSearchPtr;
typedef boost::shared_ptr<RefinementStrategy> RefinementStrategyPtr;
typedef boost::shared_ptr<RefinementSnake> RefinementSnakePtr;
typedef boost::shared_ptr<NelderMead> NelderMeadPtr;

class CSV;
class PNGFile;
class TextManager;
typedef boost::shared_ptr<PNGFile> PNGFilePtr;
typedef boost::shared_ptr<TextManager> TextManagerPtr;
typedef boost::shared_ptr<CSV> CSVPtr;

class Sampler;
typedef boost::shared_ptr<Sampler> SamplerPtr;

typedef std::vector<AtomWkr> AtomList;

class Parser;
typedef boost::shared_ptr<Parser> ParserPtr;

class VScope;
class LeftThing;
class Thing;
typedef boost::shared_ptr<VScope> VScopePtr;
typedef boost::shared_ptr<LeftThing> LeftThingPtr;
typedef boost::shared_ptr<Thing> ThingPtr;

typedef enum
{
	MaskUnchecked = 0,
	MaskProtein = 1,
	MaskEmpty = 2,

	MaskFree = 8,
	MaskWork = 9,
} MaskType;


typedef enum
{
	PDBTypeEnsemble,
	PDBTypeAverage,
	PDBTypeSamePosition,
	PDBTypeSameBFactor,
} PDBType;

typedef enum
{
	ParamOptionTorsion,
	ParamOptionKick,
	ParamOptionDampen,
	ParamOptionMagicAngles,
	ParamOptionNumBonds,
} ParamOptionType;


#endif
