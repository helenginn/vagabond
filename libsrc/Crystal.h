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

#define HARD_CODED_RESOLUTION 1.0

typedef std::map<std::string, MoleculePtr> MoleculeMap;

class Crystal : public Object, public boost::enable_shared_from_this<Crystal>, public Parser
{
public:
	Crystal();
	void addMolecule(MoleculePtr molecule);
	double concludeRefinement(int cycleNum, DiffractionPtr data);

	long int moleculeCount()
	{
		return _molecules.size();
	}

	MoleculePtr molecule(long int i)
	{
		MoleculeMap::iterator it = _molecules.begin();
		std::advance(it, i);
		return it->second;
	}

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

	mat3x3 getReal2Frac()
	{
		return _real2frac;
	}

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

	void setAnchors();
	void changeAnchors(int newAnchor);
	void tiedUpScattering();
	void realSpaceClutter(double maxRes);
	void writeMillersToFile(DiffractionPtr data, std::string prefix = "");

	void fourierTransform(int dir, double res = FLT_MAX);
	void scaleToDiffraction(DiffractionPtr data);
	void scaleComponents(DiffractionPtr data);
	void scaleSolvent(DiffractionPtr data);
	double scaleAndAddSolventScore(DiffractionPtr data = DiffractionPtr());
	double rFactorWithDiffraction(DiffractionPtr data, bool verbose = false);
	double valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
	                            bool verbose = false, double lowRes = 0,
	double highRes = 0);
	double getDataInformation(DiffractionPtr data, double partsFo = 2,
	                          double partsFc = 1, std::string prefix = "");
	void applyScaleFactor(double scale, double lowRes = 0, double highRes = 0);

	static double scaleSolventScore(void *object);
	
	static double getSolvScale(void *object)
	{
		return static_cast<Crystal *>(object)->_solvScale;
	}
	
	static void setSolvScale(void *object, double value)
	{
		static_cast<Crystal *>(object)->_solvScale = value;
	}

	static double getSolvBFac(void *object)
	{
		return static_cast<Crystal *>(object)->_solvBFac;
	}
	
	static void setSolvBFac(void *object, double value)
	{
		static_cast<Crystal *>(object)->_solvBFac = value;
	}

	void reconfigureUnitCell();
	void setupSymmetry();
	void summary();

	void makePowders();
	void tieAtomsUp();

	/* Tolerance in Angstroms. */
	std::vector<AtomPtr> getCloseAtoms(AtomPtr one, double tol = 2, bool cache = false);
	
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

	std::vector<double> getUnitCell()
	{
		return _unitCell;	
	}

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

	std::string agreementSummary();

	virtual std::string getClassName()
	{
		return "Crystal";
	}

	void fitWholeMolecules(bool translation = false, bool rotation = true);

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
	FFTPtr _solvent;
	FFTPtr _difft;
};

#endif /* defined(__vagabond__Crystal__) */
