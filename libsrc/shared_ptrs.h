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
#define SOLVENT_BITS (8 * sizeof(MaskType) / 2)

#define ATOM_SAMPLING_COUNT (30)
#define PROTEIN_SAMPLING (0.5)
#define WATER_RADIUS 0.6
#define FUTURE_MAGIC_ATOMS 12

#define deg2rad(a) ((a) * M_PI / 180)
#define rad2deg(a) ((a) / M_PI * 180)

#define b2var(a) ((a) / (8 * M_PI * M_PI))
#define var2b(a) ((a) * (8 * M_PI * M_PI))

#define ToBondPtr(a) (boost::static_pointer_cast<Bond>((a)))
#define ToBondGroupPtr(a) (boost::static_pointer_cast<BondGroup>((a)))
#define ToGhostBondPtr(a) (boost::static_pointer_cast<GhostBond>((a)))
#define ToAbsolutePtr(a) (boost::static_pointer_cast<Absolute>((a)))
#define ToAnchorPtr(a) (boost::static_pointer_cast<Anchor>((a)))
#define ToModelPtr(a) (boost::static_pointer_cast<Model>((a)))
#define ToNovalentPtr(a) (boost::static_pointer_cast<Novalent>((a)))
#define ToSpongePtr(a) (boost::static_pointer_cast<Sponge>((a)))
#define ToExplicitModelPtr(a) (boost::static_pointer_cast<ExplicitModel>((a)))
#define ToPolymerPtr(a) (boost::static_pointer_cast<Polymer>((a)))
#define ToWaterNetworkPtr(a) (boost::static_pointer_cast<WaterNetwork>((a)))
#define ToMotionPtr(a) (boost::static_pointer_cast<Motion>((a)))
#define ToMoleculePtr(a) (boost::static_pointer_cast<Molecule>((a)))
#define ToSidechainPtr(a) (boost::static_pointer_cast<Sidechain>((a)))
#define ToBackbonePtr(a) (boost::static_pointer_cast<Backbone>((a)))
#define ToCrystalPtr(a) (boost::static_pointer_cast<Crystal>((a)))
#define ToKeyGroupPtr(a) (boost::static_pointer_cast<KeyPoints>((a)))
#define ToAtomPtr(a) (boost::static_pointer_cast<Atom>((a)))
#define ToWhackPtr(a) (boost::static_pointer_cast<Whack>((a)))
#define ToTwistPtr(a) (boost::static_pointer_cast<Twist>((a)))
#define ToAtomGroupPtr(a) (boost::static_pointer_cast<AtomGroup>((a)))
#define ToMonomerPtr(a) (boost::static_pointer_cast<Monomer>((a)))
#define ToParserPtr(a) (boost::static_pointer_cast<Parser>((a)))
#define ToGridPtr(a) (boost::static_pointer_cast<RefinementGridSearch>(a))

#define ToThingPtr(a) (boost::static_pointer_cast<Thing>((a)))

class Clusterable;
class RefinableDouble;
typedef boost::shared_ptr<Clusterable> ClusterablePtr;
typedef boost::weak_ptr<Clusterable> ClusterableWkr;
typedef boost::shared_ptr<RefinableDouble> RefinableDoublePtr;

class Notifiable;
typedef boost::shared_ptr<Notifiable> NotifiablePtr;

class Options;
typedef boost::shared_ptr<Options> OptionsPtr;

class Balance;
typedef boost::shared_ptr<Balance> BalancePtr;

class FFT;
class VagFFT;
typedef boost::shared_ptr<FFT> FFTPtr;
typedef boost::shared_ptr<VagFFT> VagFFTPtr;

class Crystal;
typedef boost::shared_ptr<Crystal> CrystalPtr;
typedef boost::weak_ptr<Crystal> CrystalWkr;

class Chromosomal;
typedef boost::shared_ptr<Chromosomal> ChromosomalPtr;

class WaterNetwork;
class WaterCluster;
typedef boost::shared_ptr<WaterNetwork> WaterNetworkPtr;
typedef boost::shared_ptr<WaterCluster> WaterClusterPtr;

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

class KeyPoints;
typedef boost::shared_ptr<KeyPoints> KeyPointsPtr;
typedef boost::weak_ptr<KeyPoints> KeyPointsWkr;

class Atom;
typedef boost::shared_ptr<Atom> AtomPtr;
typedef boost::weak_ptr<Atom> AtomWkr;

class Element;
typedef boost::shared_ptr<Element> ElementPtr;

class Motion;
typedef boost::shared_ptr<Motion> MotionPtr;

class Anchor;
class Model;
class Novalent;
class Sponge;
class ExplicitModel;
class Absolute;
class Chelate;
class Bond;
class BondGroup;
class GhostBond;

typedef boost::shared_ptr<Absolute> AbsolutePtr;
typedef boost::shared_ptr<Anchor> AnchorPtr;
typedef boost::shared_ptr<Novalent> NovalentPtr;
typedef boost::shared_ptr<Sponge> SpongePtr;
typedef boost::weak_ptr<Anchor> AnchorWkr;
typedef boost::shared_ptr<Model> ModelPtr;
typedef boost::weak_ptr<Model> ModelWkr;
typedef boost::shared_ptr<Chelate> ChelatePtr;
typedef boost::shared_ptr<ExplicitModel> ExplicitModelPtr;
typedef boost::weak_ptr<ExplicitModel> ExplicitModelWkr;
typedef boost::shared_ptr<Bond> BondPtr;
typedef boost::weak_ptr<Bond> BondWkr;
typedef boost::shared_ptr<BondGroup> BondGroupPtr;
typedef boost::shared_ptr<GhostBond> GhostBondPtr;

class Whack;
class Twist;
typedef boost::shared_ptr<Whack> WhackPtr;
typedef boost::weak_ptr<Whack> WhackWkr;
typedef boost::shared_ptr<Twist> TwistPtr;
typedef boost::weak_ptr<Twist> TwistWkr;

class Bucket;
class PartialStructure;
class BucketBulkSolvent;

typedef boost::shared_ptr<Bucket> BucketPtr;
typedef boost::shared_ptr<PartialStructure> PartialStructurePtr;
typedef boost::shared_ptr<BucketBulkSolvent> BucketBulkSolventPtr;

class Dataset;
class Diffraction;
class DiffractionMtz;
typedef boost::shared_ptr<Dataset> DatasetPtr;
typedef boost::shared_ptr<Diffraction> DiffractionPtr;
typedef boost::shared_ptr<DiffractionMtz> DiffractionMtzPtr;

class Object;
typedef boost::shared_ptr<Object> ObjectPtr;

class ParamBand;
typedef boost::shared_ptr<ParamBand> ParamBandPtr;

class RefineMat3x3;
typedef boost::shared_ptr<RefineMat3x3> RefineMat3x3Ptr;

class Converter;
class RefinementGridSearch;
class RefinementStepSearch;
class RefinementStrategy;
class RefinementList;
class RefinementLBFGS;
class RefinementNelderMead;
typedef boost::shared_ptr<Converter> ConverterPtr;
typedef boost::shared_ptr<RefinementStepSearch> RefinementStepSearchPtr;
typedef boost::shared_ptr<RefinementGridSearch> RefinementGridSearchPtr;
typedef boost::shared_ptr<RefinementStrategy> RefinementStrategyPtr;
typedef boost::shared_ptr<RefinementList> RefinementListPtr;
typedef boost::shared_ptr<RefinementLBFGS> RefinementLBFGSPtr;
typedef boost::shared_ptr<RefinementNelderMead> NelderMeadPtr;

class CSV;
class PNGFile;
class TextManager;
typedef boost::shared_ptr<PNGFile> PNGFilePtr;
typedef boost::shared_ptr<TextManager> TextManagerPtr;
typedef boost::shared_ptr<CSV> CSVPtr;

class Sampler;
typedef boost::shared_ptr<Sampler> SamplerPtr;

typedef std::vector<AtomPtr> AtomList;

class BaseParser;
class Parser;
typedef boost::shared_ptr<BaseParser> BaseParserPtr;
typedef boost::weak_ptr<BaseParser> BaseParserWkr;
typedef boost::shared_ptr<Parser> ParserPtr;
typedef boost::weak_ptr<Parser> ParserWkr;

class VScope;
class LeftThing;
class Thing;
typedef boost::shared_ptr<VScope> VScopePtr;
typedef boost::shared_ptr<LeftThing> LeftThingPtr;
typedef boost::shared_ptr<Thing> ThingPtr;

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
	ParamOptionTwist,
	ParamOptionShift,
	ParamOptionBondAngle,
	ParamOptionCirclePortion,
	ParamOptionKick,
	ParamOptionOccupancy,
	ParamOptionMagicAngles,
	ParamOptionNumBonds,
	ParamOptionMaxTries,
} ParamOptionType;

typedef enum
{
	ScalingTypeShell = 0,
	ScalingTypeAbsBFactor = 1,
	ScalingTypeAbs = 2,
} ScalingType;

typedef double (*Getter)(void *);
typedef void (*Setter)(void *, double newValue);

typedef int MaskType;

#endif
