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
#include "csymlib.h"
#include <iostream>

#define HARD_CODED_RESOLUTION 1.0

typedef std::map<std::string, MoleculePtr> MoleculeMap;

class Crystal : public Object, public std::enable_shared_from_this<Crystal>
{
public:
	Crystal();
	void addMolecule(MoleculePtr molecule);
	void concludeRefinement(int cycleNum, DiffractionPtr data,
							CrystalPtr crystal);

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
	void realSpaceClutter();
	void writeCalcMillersToFile(DiffractionPtr data, std::string prefix = "");

	void fourierTransform(int dir);
	void scaleToDiffraction(DiffractionPtr data);
	double rFactorWithDiffraction(DiffractionPtr data, bool verbose = false);
	double valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
								bool verbose = false, double lowRes = 0,
								double highRes = 0);
	void getDataInformation(DiffractionPtr data, double partsFo = 2,
							  double partsFc = 1);
	void applyScaleFactor(double scale, double lowRes = 0, double highRes = 0);

	void reconfigureUnitCell();
	void summary();

	void tieAtomsUp();
	
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

	void setOverallBFactor(double b)
	{
		_overallB = b;
	}

	double getOverallBFactor()
	{
		return _overallB;
	}

private:
	MoleculeMap _molecules;
	std::string _filename;

	std::vector<double> _unitCell;
	double _firstScale;
	mat3x3 _hkl2real;
	mat3x3 _real2frac;
	CSym::CCP4SPG *_spaceGroup;
	double _maxResolution;
	std::vector<int> _anchorResidues;
	double _overallB;
	double totalToScale();

	void applySymOps();

	FFTPtr _fft;
	FFTPtr _difft;
};

#endif /* defined(__vagabond__Crystal__) */
