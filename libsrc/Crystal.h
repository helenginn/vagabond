// Vagabond : bond-based macromolecular model refinement
// Copyright (C) 2017-2018 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.
//

#ifndef __vagabond__Crystal__
#define __vagabond__Crystal__

#include <stdio.h>
#include <vector>
#include "shared_ptrs.h"
#include <hcsrc/mat3x3.h>
#include "Object.h"
#include <string>
#include <map>
#include <hcsrc/maths.h>
#include "Molecule.h"
#include "../libccp4/csymlib.h"
#include <iostream>

/**
 * \class Crystal
 * \brief The concept of a crystal arranged in a certain lattice containing
 * various molecules and solvent.
 *
 * Crystal looks after the calculation of Fcs and also compares against a
 * separate DiffractionData object, but it endeavours to be separate from the
 * DiffractionData conceptually.
 */

typedef enum
{
	JobPositions,
	JobWholeMol,
	JobIntraMol,
	JobBackbone,
	JobSidechain,
} JobType;

typedef std::map<std::string, MoleculePtr> MoleculeMap;

typedef double (*option_getter)();

typedef struct
{
	double minRes;
	double maxRes;
	double rFactor;
	double scale;
	double std_err;
	double aveFo;
	double phi_spread;
	double count;
	std::vector<double> work1;
	std::vector<double> work2;
	std::vector<double> free1;
	std::vector<double> free2;
} ShellInfo;

inline ShellInfo makeShellInfo(double min, double max)
{
	ShellInfo shell;
	shell.minRes = min;
	shell.maxRes = max;
	shell.rFactor = 0;
	shell.scale = 1;
	shell.std_err = 0;
	shell.aveFo = 0;
	shell.phi_spread = 0;
	shell.count = 0;
	
	return shell;
}

/* resolutions are in real space */
inline int findShell(std::vector<ShellInfo> &shells, double res)
{
	if (res < shells.back().maxRes || 
	    res > shells.front().minRes)
	{
		return -1;
	}
	
	int size = shells.size();
	int min = 0;
	int max = size - 1;
	
	if (res >= shells[min].maxRes) return min;
	if (res <= shells[max].minRes) return max;
	
	while (true)
	{
		int chop = (max + min) / 2;
		int higher = (res <= shells[chop].maxRes);
		if (higher) 
		{
			min = chop;
		}
		else 
		{
			max = chop;
		}
		
		if (max - min == 1)
		{
			return max;
		}
	}
}

inline void makeShells(std::vector<ShellInfo> *shells, double min, double max,
                       int number = 20)
{
	/* Then apply to individual resolution bins */
	std::vector<double> bins;
	generateResolutionBins(min, max, number, &bins);

	/* Extend the final bin by a little bit, so as not to lose any
	 * stragglers. */
	bins[bins.size() - 1] *= 0.95;

	/* Make the series of shells */
	shells->clear();

	for (size_t i = 0; i < bins.size() - 1; i++)
	{
		ShellInfo shell = makeShellInfo(bins[i], bins[i + 1]);
		shells->push_back(shell);
	}

}


namespace Vagabond
{
class Crystal : public Object, public AtomGroup
{
public:
	CrystalPtr shared_from_this()
	{
		return ToCrystalPtr(BaseParser::shared_from_this());
	}
	
	Crystal();
	virtual ~Crystal() {};
	void addMolecule(MoleculePtr molecule);
	
	/**
	* 	Conclude refinement finishes a round of refinement and triggesr an
	* 	FFT, PDB output, graphing and statistics.
	* 	\param cycleNum cycle number to output filenames appropriately
	* 	\param data diffraction data against which statistics should be generated.
	*/
//	double concludeRefinement(int cycleNum, DiffractionPtr data);
	
	/** Should be folded into previous concludeRefinement(...) soon */
	void wrapUpRefinement();
	
	/** Move back to an earlier saved version of the model */
	static void vsRestoreState(void *object, double val);

	void scaleToDiffraction(DiffractionPtr data, bool full = true);

	bool isSilent()
	{
		return _silent;
	}
	
	void setSilent(bool silent)
	{
		_silent = silent;
	}
	
	/**
	* How many molecules are included in a Crystal.
	*/
	size_t moleculeCount()
	{
		return _molecules.size();
	}

	/**
	* 	get the stored molecule. Cannot guarantee order of molecules will be
	* 	the same as in the PDB file, but adequate for looping over molecules.
	*   \param i molecule number to obtain.
	* 	*/
	MoleculePtr molecule(long int i)
	{
		MoleculeMap::iterator it = _molecules.begin();
		std::advance(it, i);
		return it->second;
	}

	/**
	* 	get the stored molecule by chain name.
	*/
	MoleculePtr molecule(std::string chain)
	{
		if (_molecules.count(chain))
		{
			return _molecules[chain];
		}

		return MoleculePtr();
	}

	void recalculateAtoms();

	void setReal2Frac(mat3x3 mat);

	void setHKL2Real(mat3x3 mat);

	/**
	* 	Returns matrix turning real coordinates (of atoms) to fractional
	* 	coordinates of the unit cell. */
	mat3x3 getReal2Frac()
	{
		return _real2frac;
	}

	/**
	* 	Returns matrix turning Miller index coordinates to fractional
	* 	coordinates of the reciprocal cell. */
	mat3x3 getFrac2Real()
	{
		return _hkl2real;
	}

	VagFFTPtr getFFT()
	{
		return _fft;
	}
	
	void clearFFT()
	{
		_fft = VagFFTPtr();
		_sampling = -1;
	}

	VagFFTPtr getDiFFT()
	{
		return _difft;
	}

	VagFFTPtr getOrigDensity()
	{
		return _original;
	}
	
	size_t symOpCount()
	{
		return _spaceGroup->nsymop;
	}
	
	void prepareFFT(VagFFTPtr ft);
	/** Write out the % contribution of elements in the crystal to the
	 *  total difference density */
	void differenceAttribution();
	
	void removeAtom(AtomPtr atom);
	void removeMolecule(MoleculePtr mol);

	/** Calculates the anchor residue for each Polymer and assigns to each. */
	void setAnchors();
	
	size_t polymerCount();
	
	/** Prints scattering proportion in this crystal determined by bonds. */
	void tiedUpScattering();
	
	double averageBFactor();
	void scaleAnchorBs(double ratio);
	
	/** Loops through all molecules and places them in the map.
	* 	\param maxRes max resolution used to determine voxel spacing if it has
	* 	not yet been determined, in Angstroms. */
	void realSpaceClutter();
	
	double rFactorWithDiffraction(DiffractionPtr data, bool verbose = false);
	double valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
	                            bool allShells = false, bool verbose = false);
	
	void scaleAndBFactor(DiffractionPtr data, double *scale, 
                              double *bFactor, VagFFTPtr model = VagFFTPtr());

	/** Applies scale factor within the resolution ranges.
	* 	\param scale factor (multiplies Fc with this)
	* 	\param lowRes low resolution in Angstroms, 0 ignores this
	* 	\param highRes high resolution in Angstroms, 0 ignores this
	*/
	void applyScaleFactor(double scale, double lowRes = 0, double highRes = 0,
	                      double bFactor = 0);

	void applyShellFactors(DiffractionPtr data);
	double getAdjustBFactor();
	
	void undo();
	bool undoIfWorse();
	bool returnToBestState();

	void setupSymmetry();
	void summary();

	void tieAtomsUp();
	
	void setLastLocalCC(double cc)
	{
		_lastLocalCC = _localCC;
		_localCC = cc;
		if (_bestLocalCC < cc)
		{
			_bestLocalCC = cc;
		}
	}
	
	void hydrogenateContents();
	
	/**
	* 	
	*/
	AtomPtr getClosestAtom(vec3 pos);

	void setFilename(std::string file)
	{
		_filename = file;
	}

	std::string getFilename()
	{
		return _filename;
	}

	void setSpaceGroup(CSym::CCP4SPG *spg)
	{
		_spaceGroup = spg;
	}

	CSym::CCP4SPG *getSpaceGroup()
	{
		return _spaceGroup;
	}
	
	/**	Set unit cell dimensions of Crystal. Nothing else is updated. */
	void setUnitCell(double a, double b, double c,
	                 double alpha, double beta, double gamma)
	{
		_unitCell.clear();
		_unitCell.push_back(a);
		_unitCell.push_back(b);
		_unitCell.push_back(c);
		_unitCell.push_back(alpha);
		_unitCell.push_back(beta);
		_unitCell.push_back(gamma);
	}

	/** Returns unit cell for Crystal
 * 	\return vector with 6 params, three dimensions in Angstroms, three angles
 * 	in degrees. */
	std::vector<double> getUnitCell()
	{
		return _unitCell;	
	}

	/** Set maximum resolution
	* 	\param maxRes resolution in Angstroms.
	* 	*/
	void setMaxResolution(double maxRes)
	{
		_maxResolution = maxRes;
	}

	int totalAnchors()
	{
		return _anchorResidues.size();
	}

	/** Returns a string describing the R factors and CCs.
	* 	\return summary string in human-readable format. */
	std::string agreementSummary();
	
	virtual std::string getClassName()
	{
		return "Crystal";
	}
	
	void makeBucket();
	
	BucketPtr getBucket()
	{
		return _bucket;
	}
	
	double getDataWilsonB()
	{
		return _dataWilsonB;
	}

	bool applyWilsonToAnchors();

	void rigidBodyRefinement();
	/** Fit whole molecules to electron density as a refinement protocol.
	* 	\param translation if true, will refine translation parameters
	* 	\param rotation if true, will refine rotation parameters. */
	void fitWholeMolecules(bool recip = false);
	void refinePolymers(RefinementType type);
	void refinePositions(int total = 15);
	bool refineThreaded(JobType type, int total = 15);
	void refineSidechainPositions();
	void refineSidechains();
	void refineCrude();
	void pruneWaters();
	void resetMotions();
	void savePositions();
	void refreshAnchors();
	void closenessSummary();
	void makeOverallMotion();
	MotionPtr getOverallMotion();
	void reindexAtoms(mat3x3 reindex, vec3 trans);

	void addPDBContents(std::string pdb);
	void updatePDBContents(std::string pdbName);
	
	void addMotion(MotionPtr mot, PolymerPtr pol = PolymerPtr());
	
	size_t motionCount()
	{
		return _motions.size();
	}
	
	double getRWork()
	{
		return _rWork;
	}
	
	double getCCWork()
	{
		return _ccWork;
	}
	
	double getWorkValue()
	{
		return _rWork;
	}
	
	int getCycleNum()
	{
		return _cycleNum;
	}
	
	void addComment(std::string comment)
	{
		_comments += "(" + _vbondFile + ") ";
		_comments += comment;
		_comments += "\n";
	}
	
	WaterNetworkPtr getWaterNetwork();
	
	void openInCoot();
	
	/** Obtain the current number of samples, i.e. number of conformers
	 * in the ensemble. */
	int getSampleNum();
	
	/** Obtain the current B factor applied to any many-positioned models of
	 * atoms. */
	double getRealBFactor();
	
	double getProbeRadius();

	std::vector<AtomPtr> getHydrogenBonders();
	
	/** Returns the maximum resolution. Diffraction data as the input
	 * target is required; resolution will be determined from this if not
	 * already determined, or input as command line option. */
	double getMaxResolution(DiffractionPtr data = DiffractionPtr());

	double getProteinSampling();
	
	void addBlob(BlobPtr blob);
	void removeBlob(BlobPtr blob);
	
	void updateLargestNum(AtomPtr atom);
	void whack();
	void chelate();

	int issueAtomNumber()
	{
		if (_largestNum == -INT_MAX)
		{
			return 0;
		}
		return _largestNum + 1;
	}
	
	void setSampleNum(int num)
	{
		_sampleNum = num;
	}
	
	std::vector<ShellInfo> getShells()
	{
		return _shells;
	}
	
	/** Scale any data set that has been provided as FPART/PHIPART in
	 *  the input file */
	void scaleAnyPartialSet();

	void bestGlobalParameters();
	void writeVagabondFile(int cycleNum);
	void makePDBs(std::string suffix);

	void applySymOps();
	void fourierTransform(int dir, double res = -1);
protected:
	virtual void postRestoreState();
	virtual void addObject(ParserPtr object, std::string category);
	virtual std::string getParserIdentifier()
	{
		return "Crystal_" + _filename;
	}
	
	virtual void addProperties();
	virtual void postParseTidy();
private:
	void addMissingAtoms(std::vector<AtomPtr> atoms);
	void saveAtomsForThreading();
	void fusePolymers();

	MoleculeMap _molecules;
	std::string _filename;
	std::string _vbondFile;

	std::vector<double> _unitCell;
	mat3x3 _hkl2real;
	mat3x3 _real2frac;
	CSym::CCP4SPG *_spaceGroup;
	int _spgNum;
	std::string _spgString;
	bool _tied;

	double _maxResolution;
	std::vector<int> _anchorResidues;
	double totalToScale();
	double _calcElec;

	void reportScaling();
	void setupOriginalMap();
	
	std::string _lastEnsemblePDB;
	std::string _lastAveragePDB;
	std::string _lastMtz;
	std::string _comments;
	
	std::map<double, double> _resBinAves;
	std::vector<ShellInfo> _shells;
	std::vector<MotionPtr> _motions;

	double _rWork;
	double _rFree;
	double _ccWork;
	double _ccFree;
	double _lastMetric;
	double _bestMetric;
	double _localCC;
	double _lastLocalCC;
	double _bestLocalCC;
	
	double _realBFactor;
	double _probeRadius;
	double _sampling;

	int _sinceBestNum;
	int _bestState;
	int _cycleNum;
	int _sampleNum;
	bool _silent;
	
	double _bFacFit;
	int _correlPlotNum;
	DiffractionPtr _data;

	VagFFTPtr _fft;
	VagFFTPtr _difft;
	
	/* imag component may contain (weighted map - original) */
	VagFFTPtr _original;

	BucketPtr _bucket;
	int _largestNum;
	double _dataWilsonB;

	double updateVariable(double *local, option_getter get, Setter set,
	                      std::string name, std::string unit, 
	                      double default_val);
};
}

#endif /* defined(__vagabond__Crystal__) */
