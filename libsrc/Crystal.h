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

typedef std::map<std::string, MoleculePtr> MoleculeMap;

class Crystal : public Object
{
public:
	void addMolecule(MoleculePtr molecule);

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

	void setReal2HKL(mat3x3 mat);
	void setHKL2Real(mat3x3 mat);
	mat3x3 getReal2Frac()
	{
		return _real2frac;
	}

	FFTPtr getFFT()
	{
		return _fft;
	}

	void realSpaceClutter();
	void writeCalcMillersToFile(double resolution = 1.0);

	void fourierTransform(int dir)
	{
		_fft->fft(dir);
	}

	void scaleToDiffraction(DiffractionPtr data);
	double rFactorWithDiffraction(DiffractionPtr data, bool verbose = false);
	double valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
								bool verbose = false);
	void transplantAmplitudes(DiffractionPtr data, double partsFo = 2,
							  double partsFc = 1);
	
	void summary();

	void tieAtomsUp();
	
	void setFilename(std::string file)
	{
		_filename = file;
	}
	
private:
	MoleculeMap _molecules;
	std::string _filename;


	mat3x3 _hkl2real;
	mat3x3 _real2frac;

	FFTPtr _fft;
};

#endif /* defined(__vagabond__Crystal__) */
