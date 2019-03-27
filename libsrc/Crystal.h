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
#include "mat3x3.h"
#include "Object.h"
#include "fftw3d.h"
#include <string>
#include <map>
#include "maths.h"
#include "Molecule.h"
#include "../libccp4/csymlib.h"
#include <iostream>
#include "Parser.h"

/**
 * \class Crystal
 * \brief The concept of a crystal arranged in a certain lattice containing
 * various molecules and solvent.
 *
 * Crystal looks after the calculation of Fcs and also compares against a
 * separate DiffractionData object, but it endeavours to be separate from the
 * DiffractionData conceptually.
 */

typedef std::map<std::string, MoleculePtr> MoleculeMap;

typedef struct
{
	double minRes;
	double maxRes;
	double rFactor;
	double scale;
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
	
	return shell;
}

inline double getNLimit(FFTPtr fftData, FFTPtr fftModel, int axis = 0)
{
	double nLimit = std::min(*(&fftData->nx + axis), 
	                         *(&fftModel->nx + axis));

	nLimit = nLimit - ((int)nLimit % 2);
	nLimit /= 2;

	return nLimit;	
}

inline vec3 getNLimits(FFTPtr data, FFTPtr fftModel)
{
	vec3 lims;
	lims.x = getNLimit(data, fftModel, 0);
	lims.y = getNLimit(data, fftModel, 1);
	lims.z = getNLimit(data, fftModel, 2);
	return lims;
}

class Crystal : public Object, public Parser
{
public:
	CrystalPtr shared_from_this()
	{
		return ToCrystalPtr(Parser::shared_from_this());
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
	double concludeRefinement(int cycleNum, DiffractionPtr data);
	
	/** Should be folded into previous concludeRefinement(...) soon */
	static double vsConcludeRefinement(void *object);
	
	/** Calculate new observed/calculated density but don't write out
	 * 	R factors to the screen */
	void silentConcludeRefinement();
	
	/** Move back to an earlier saved version of the model */
	static void vsRestoreState(void *object, double val);

	bool isSilent()
	{
		return _silent;
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
	mat3x3 getHKL2Frac()
	{
		return _hkl2real;
	}

	FFTPtr getFFT()
	{
		return _fft;
	}

	FFTPtr getDiFFT()
	{
		return _difft;
	}
	
	size_t symOpCount()
	{
		return _spaceGroup->nsymop;
	}
	
	void omitScan();
	
	/** Write out the % contribution of elements in the crystal to the
	 *  total difference density */
	void differenceAttribution();
	
	void removeAtom(AtomPtr atom);
	
	/** Takes an input crystal position and moves it to the nearest
	 *  position which falls on an integer value of the FFT grid */
	vec3 snapToGrid(vec3 pos);

	/** Calculates the anchor residue for each Polymer and assigns to each. */
	void setAnchors();
	
	/** Prints scattering proportion in this crystal determined by bonds. */
	void tiedUpScattering();
	
	double averageBFactor();
	void scaleAnchorBs(double ratio);
	
	/** Loops through all molecules and places them in the map.
	* 	\param maxRes max resolution used to determine voxel spacing if it has
	* 	not yet been determined, in Angstroms. */
	void realSpaceClutter(double maxRes);
	
	/** Creates an MTZ file to open in Coot.
	* 	\param data diffraction data to get F-obs from.
	* 	\param prefix for choosing a filename */
	void writeMillersToFile(DiffractionPtr data, std::string prefix = "");

	void fourierTransform(int dir, double res = FLT_MAX);
	
	/** Finds appropriate scale factor for the solvent contribution by means
	 * of a grid search. Starts with a global scale factor, then scales
	 * the solvent, then a full scale by resolution bin. */
	void scaleComponents(DiffractionPtr data);
	double rFactorWithDiffraction(DiffractionPtr data, bool verbose = false);
	double valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
	                            bool allShells = false, bool verbose = false);
	double getDataInformation(DiffractionPtr data, double partsFo = 2,
	                          double partsFc = 1, std::string prefix = "");
	
	void scaleAndBFactor(DiffractionPtr data, double *scale, 
                              double *bFactor, FFTPtr model = FFTPtr());

	/** Applies scale factor within the resolution ranges.
	* 	\param scale factor (multiplies Fc with this)
	* 	\param lowRes low resolution in Angstroms, 0 ignores this
	* 	\param highRes high resolution in Angstroms, 0 ignores this
	*/
	void applyScaleFactor(double scale, double lowRes = 0, double highRes = 0,
	                      double bFactor = 0);

	double applyShellFactors(DiffractionPtr data);
	double getAdjustBFactor();
	
	bool undoIfWorse();

	void setupSymmetry();
	void summary();

	void tieAtomsUp();
	
	void hydrogenateContents();
	
	/**
	* Find all close atoms within this crystal to a chosen atom. Including
	* the chosen atom.
	* \param one chosen atom.
	* \param tol tolerance in Angstroms - furthest point from a given atom.
	* \param cache if true, instead of returning atoms, cache given atoms for
	* quicker checking routine. Until cache is cleared, this subset will be
	* searched in the future.
	*/
	std::vector<AtomPtr> getCloseAtoms(AtomPtr one, double tol = 2, 
	                                   bool cache = false);
	
	std::vector<AtomPtr> getAtomsInBox(vec3 target, double tolx,
	                                   double toly, double tolz);

	
	std::vector<AtomPtr> getCloseAtoms(std::vector<AtomPtr> atoms,
	                                   double tol = 2); 

	/**
	* 	
	*/
	AtomPtr getClosestAtom(vec3 pos);
	
	/**
	* 	Clear close atom cache and research the entire Crystal next time
	* 	getCloseAtoms is called. */
	void clearCloseCache();

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
	
	BucketPtr getBucket()
	{
		return _bucket;
	}

	/** Fit whole molecules to electron density as a refinement protocol.
	* 	\param translation if true, will refine translation parameters
	* 	\param rotation if true, will refine rotation parameters. */
	void fitWholeMolecules();
	void refinePolymers(RefinementType type);
	void refinePositions();
	bool refineIntraMovements();
	void refineSidechainPositions();
	void refineSidechains();
	void refineCrude();
	void savePositions();
	void rigidBodyFit();
	
	double getRWork()
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
	
	double getMaximumDStar(DiffractionPtr data = DiffractionPtr());
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
	double getMaxResolution(DiffractionPtr data);
	
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
	
	/** Scale any data set that has been provided as FPART/PHIPART in
	 *  the input file */
	void scaleAnyPartialSet();
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

	void makePDBs(std::string suffix);
	void writeVagabondFile(int cycleNum);
	void applySymOps();
	
	std::string _lastEnsemblePDB;
	std::string _lastAveragePDB;
	std::string _lastMtz;
	std::string _comments;
	
	std::map<double, double> _resBinAves;
	std::vector<ShellInfo> _shells;

	double _rWork;
	double _rFree;
	double _ccWork;
	double _ccFree;
	double _lastRWork;
	double _bestRWork;
	
	double _realBFactor;
	double _probeRadius;

	int _sinceBestNum;
	int _bestState;
	int _cycleNum;
	int _sampleNum;
	bool _silent;
	
	double _solvScale;
	double _solvBFac;
	double _bFacFit;
	int _correlPlotNum;
	DiffractionPtr _data;

	FFTPtr _fft;
	FFTPtr _difft;

	BucketPtr _bucket;
	int _largestNum;
	
	void scaleSolvent(DiffractionPtr data);
	void scaleToDiffraction(DiffractionPtr data, bool full = true);
};

#endif /* defined(__vagabond__Crystal__) */
