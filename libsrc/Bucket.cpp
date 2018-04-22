//
//  Bucket.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bucket.h"
#include "Atom.h"
#include "Node.h"
#include <iostream>
#include <iomanip>
#include "fftw3d.h"
#include "RefinementGridSearch.h"
#include "Shouter.h"
#include "Crystal.h"
#include "Diffraction.h"
#include "CSV.h"
#include "Options.h"

#define MAX_CHECK_DISTANCE 6.0
#define MIN_CHECK_DISTANCE 1.0
#define CHECK_DISTANCE_STEP 0.30

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
	
	/** If we are doing powder analysis we don't actually want
	* 	to add the solvent */
	if (Options::shouldPowder())
	{
		return;
	}	
	
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

		_maskedRegions = FFTPtr(new FFT(*_solvent));
		processMaskedRegions();
		_solvent->fft(dir);
		applySymOps(spg, res);
		_solvent->normalise();
	}
}

void Bucket::processMaskedRegions()
{
	_maskedRegions->setupMask();
	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	mat3x3 frac2Real = mat3x3_inverse(real2Frac);

	for (long i = 0; i < _maskedRegions->nn; i++)
	{
		if (_maskedRegions->data[i][0] > 0.5)
		{
			_maskedRegions->setMask(i, 1);
			continue;
		}

		vec3 newfrac = _maskedRegions->fracFromElement(i);
		vec3 pos = mat3x3_mult_vec(frac2Real, newfrac);
		
		AtomPtr atom = getCrystal()->getClosestAtom(pos);
		
		if (atom->isHeteroAtom())
		{
			_maskedRegions->setMask(i, 2);
		}
	}
}

void Bucket::abandonCalculations()
{
	_solvent = FFTPtr();
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

/* Left, centre and right are the vec3s of interest */
void Bucket::populateHistogram(Node *node, vec3 centre, vec3 left)
{
	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	vec3 ldiff = vec3_subtract_vec3(left, centre);
	double step = CHECK_DISTANCE_STEP;

	FFTPtr fft = getCrystal()->getFFT();

	vec3 transCentre = mat3x3_mult_vec(real2Frac, centre);
	vec3 transLeft = mat3x3_mult_vec(real2Frac, left);
	double centreDensity = fft->getRealFromFrac(transCentre);
	double leftDensity = fft->getRealFromFrac(transLeft);
	
	for (double x = -MAX_CHECK_DISTANCE;
	     x <= MAX_CHECK_DISTANCE; x += step)
	{
		for (double y = -MAX_CHECK_DISTANCE;
		     y <= MAX_CHECK_DISTANCE; y += step)
		{
			for (double z = -MAX_CHECK_DISTANCE;
			     z <= MAX_CHECK_DISTANCE; z += step)
			{
				vec3 offset = make_vec3(x, y, z);
				double rlength = vec3_length(offset);
				if (rlength < MIN_CHECK_DISTANCE ||
				    rlength > MAX_CHECK_DISTANCE)
				{
					continue;
				}

				vec3 right = vec3_add_vec3(centre, offset);
				
				double angle = vec3_angle_with_vec3(ldiff, offset);
				double degrees = rad2deg(angle);
				
				if (degrees < 0) degrees = - degrees;
				
				vec3 transRight = mat3x3_mult_vec(real2Frac, right);
				double rightDensity = fft->getRealFromFrac(transRight);

				double mult = (leftDensity * rightDensity) * centreDensity;
				
				if (mult != mult || degrees != degrees || rlength != rlength)
				{
					continue;	
				}
				
				add_to_node(node, rlength, degrees, mult);
			}
		}
	}

}

void Bucket::addAnalysisForSolventPos(Node *node, vec3 centre, double distance)
{
	double step = CHECK_DISTANCE_STEP;
	double minDist = distance - step / 3;
	double maxDist = distance + step / 3;

	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	vec3 transCentre = mat3x3_mult_vec(real2Frac, centre);

	double val = _maskedRegions->getRealFromFrac(transCentre);

	if (val <= 0.8)
	{
		return;
	}
	
	double centreDensity = getCrystal()->getFFT()->getRealFromFrac(transCentre);
	
	if (centreDensity < 0)
	{
		return;
	}
	
	std::cout << "Adding data from " << vec3_desc(centre) << std::endl;
	for (double x = -maxDist - step; x <= maxDist + step; x += step)
	{
		for (double y = -maxDist - step; y <= maxDist + step; y += step)
		{
			for (double z = -maxDist - step; z <= maxDist + step; z += step)
			{
				vec3 offset = make_vec3(x, y, z);
				vec3 left = vec3_add_vec3(centre, offset);
				
				/* Temporary removal of some distances */
	
				vec3 ldiff = vec3_subtract_vec3(left, centre);
				double llength = vec3_length(ldiff);
				if (llength < minDist || llength > maxDist)
				{
					continue;
				}

				populateHistogram(node, centre, left);
			}
		}
	}

	CSVPtr csv = CSV::nodeToCSV(node);
	csv->writeToFile("solvent_density_analysis.csv");
}

void Bucket::analyseSolvent(double distance)
{
	Node *node = malloc_node();
	prepare_node(node, 5, 0, 0, MAX_CHECK_DISTANCE + 1, 180);

	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	mat3x3 frac2Real = mat3x3_inverse(real2Frac);

	for (int i = 0; i < _maskedRegions->nn; i+=2)
	{
		if (_maskedRegions->data[i][0] < 0.8) continue;

		vec3 newfrac = getCrystal()->getFFT()->fracFromElement(i);
		vec3 centre = mat3x3_mult_vec(frac2Real, newfrac);

		addAnalysisForSolventPos(node, centre, distance);
	}

	free_node(node);
}



