//
//  Crystal.h
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
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
 * Crystal looks after the calculation of Fcs and also compares against a
 * separate DiffractionData object, but it endeavours to be separate from the
 * DiffractionData conceptually.
 */

typedef std::map<std::string, MoleculePtr> MoleculeMap;

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

	/**
	* How many molecules are included in a Crystal.
	*/
	long int moleculeCount()
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
	mat3x3 getHKL2Real()
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
	
	
	void multiplyMap(double scale);

	/** Calculates the anchor residue for each Polymer and assigns to each. */
	void setAnchors();
	
	/** Prints scattering proportion in this crystal determined by bonds. */
	void tiedUpScattering();
	
	/** Loops through all molecules and places them in the map.
	* 	\param maxRes max resolution used to determine voxel spacing if it has
	* 	not yet been determined, in Angstroms. */
	void realSpaceClutter(double maxRes);
	
	/** Creates an MTZ file to open in Coot.
	* 	\param data diffraction data to get F-obs from.
	* 	\param prefix for choosing a filename */
	void writeMillersToFile(DiffractionPtr data, std::string prefix = "");

	void fourierTransform(int dir, double res = FLT_MAX);
	void scaleComponents(DiffractionPtr data);
	double scaleAndAddSolventScore(DiffractionPtr data = DiffractionPtr());
	double rFactorWithDiffraction(DiffractionPtr data, bool verbose = false);
	double valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
	                            bool verbose = false, double lowRes = 0,
	double highRes = 0);
	double getDataInformation(DiffractionPtr data, double partsFo = 2,
	                          double partsFc = 1, std::string prefix = "");
	
	/** Applies scale factor within the resolution ranges.
	* 	\param scale factor (multiplies Fc with this)
	* 	\param lowRes low resolution in Angstroms, 0 ignores this
	* 	\param highRes high resolution in Angstroms, 0 ignores this
	*/
	void applyScaleFactor(double scale, double lowRes = 0, double highRes = 0);

	void reconfigureUnitCell();
	void setupSymmetry();
	void summary();

	void makePowders();
	void tieAtomsUp();
	
	void hydrogenateContents();
	
	/**
	* Find all close atoms within this crystal to a chosen atom.
	* \param one chosen atom.
	* \param tol tolerance in Angstroms - furthest point from a given atom.
	* \param cache if true, instead of returning atoms, cache given atoms for
	* quicker checking routine. Until cache is cleared, this subset will be
	* searched in the future.
	*/
	std::vector<AtomPtr> getCloseAtoms(AtomPtr one, double tol = 2, bool cache = false);
	
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

	void addAnchorResidue(int anchor)
	{
		_anchorResidues.push_back(anchor);
		std::cout << "Adding anchor residue " << anchor << " to "
		<< getFilename() << "." << std::endl;

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

	/** Fit whole molecules to electron density as a refinement protocol.
	* 	\param translation if true, will refine translation parameters
	* 	\param rotation if true, will refine rotation parameters. */
	void fitWholeMolecules(bool translation = false, bool rotation = true);
	
	
	/** Analysis of backbone density and generation of heuristics */
	void backboneDensityAnalysis();
protected:
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

	std::vector<double> _unitCell;
	mat3x3 _hkl2real;
	mat3x3 _real2frac;
	CSym::CCP4SPG *_spaceGroup;
	int _spgNum;
	bool _tied;

	double _maxResolution;
	std::vector<int> _anchorResidues;
	double totalToScale();

	void makePDBs(std::string suffix);
	void writeVagabondFile(int cycleNum);
	void applySymOps(double res = FLT_MAX);

	double _rWork;
	double _rFree;
	double _ccWork;
	double _ccFree;
	
	double _solvScale;
	double _solvBFac;
	int _correlPlotNum;
	DiffractionPtr _data;

	FFTPtr _fft;
	FFTPtr _difft;

	BucketPtr _bucket;
	
	void scaleSolvent(DiffractionPtr data);
	void scaleToDiffraction(DiffractionPtr data);
};

#endif /* defined(__vagabond__Crystal__) */
