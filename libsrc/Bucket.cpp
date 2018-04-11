//
//  Bucket.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bucket.h"
#include <iostream>
#include "fftw3d.h"
#include "RefinementGridSearch.h"
#include "Shouter.h"
#include "Crystal.h"
#include "Diffraction.h"

void Bucket::scaleSolvent()
{
	if (!_data)
	{
		shout_at_helen("Need diffraction data to scale solvent");
	}
	
	setSolvScale(this, 4.0);
	setSolvBFac(this, 400);
	
	RefinementStrategyPtr grid;
	grid = RefinementStrategyPtr(new RefinementGridSearch());
	grid->setJobName("solvent_scale_grid_search");
	grid->setEvaluationFunction(scaleSolventScore, this);
	grid->addParameter(this, getSolvScale, setSolvScale, 8.0, 0.4, "scale");
	grid->addParameter(this, getSolvBFac, setSolvBFac, 800, 40.0, "bfac");
	grid->refine();
	
	/** Now add this into the FFT */
	
	FFTPtr fft = getCrystal()->getFFT();
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	double nLimit = fft->nx;

	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = -nLimit; k < nLimit; k++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-2 * _solvBFac / four_d_sq);

				long nModel = fft->element(i, j, k);
				float realProtein = fft->data[nModel][0];
				float imagProtein = fft->data[nModel][1];
				float realSolvent = _solvent->data[nModel][0];	
				float imagSolvent = _solvent->data[nModel][1];	
				realSolvent *= _solvScale * bFacMod;
				imagSolvent *= _solvScale * bFacMod;

				float real = realProtein + realSolvent;
				float imag = imagProtein + imagSolvent;

				fft->data[nModel][0] = real;
				fft->data[nModel][1] = imag;
			}
		}
	}
}

double Bucket::scaleSolventScore(void *object)
{
	return static_cast<Bucket *>(object)->scaleAndAddSolventScore();
}

double Bucket::scaleAndAddSolventScore()
{
	FFTPtr fftData = _data->getFFT();
	FFTPtr fft = getCrystal()->getFFT();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	double nLimit = std::min(fftData->nx, fft->nx);
	nLimit /= 2;
	mat3x3 real2frac = getCrystal()->getReal2Frac();

	std::vector<double> fData, fModel;

	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = -nLimit; k < nLimit; k++)
			{
				bool isRfree = (fftData->getMask(i, j, k) == 0);

				if (isRfree) continue;

				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(real2frac, &ijk);
				double length = vec3_length(ijk);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-2 * _solvBFac / four_d_sq);

				int _i = 0; int _j = 0; int _k = 0;	
				CSym::ccp4spg_put_in_asu(spg, i, j, k,
				                         &_i, &_j, &_k);
				float ref = sqrt(fftData->getIntensity(_i, _j, _k));

				if (ref != ref) continue;

				long nModel = fft->element(i, j, k);
				float realProtein = fft->data[nModel][0];
				float imagProtein = fft->data[nModel][1];

				float realSolvent = _solvent->data[nModel][0];	
				float imagSolvent = _solvent->data[nModel][1];	
				realSolvent *= _solvScale * bFacMod;
				imagSolvent *= _solvScale * bFacMod;

				float real = realProtein + realSolvent;
				float imag = imagProtein + imagSolvent;
				float amp = sqrt(real * real + imag * imag);

				fData.push_back(ref);
				fModel.push_back(amp);
			}
		}
	}

	double correl = correlation(fData, fModel);
	
	return -correl;
}


void Bucket::applySymOps(CSym::CCP4SPG *spaceGroup, double res)
{
	if (spaceGroup->spg_num == 1)
	{
		return;
	}

	_solvent->applySymmetry(spaceGroup, res);
}

void Bucket::fourierTransform(int dir, double res)
{
	/* Only care about reciprocal space */
	if (dir == 1)
	{
		CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();

		_solvent->fft(dir);
		applySymOps(spg, res);
		_solvent->normalise();
	}

}

void Bucket::writeMillersToFile(std::string prefix, double maxRes)
{
	std::string solventFileOnly = prefix + "_solvent_vbond.mtz";
	CrystalPtr crystal = getCrystal();
	std::vector<double> unitCell = crystal->getUnitCell();
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	
	_solvent->writeReciprocalToFile(solventFileOnly, maxRes,
	                                spg, unitCell, real2frac);
	
}
