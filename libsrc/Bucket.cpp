//
//  Bucket.cpp
//  vagabond
//
//  Created by Helen Ginn on 17/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Bucket.h"
#include "FileReader.h"
#include "Atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "fftw3d.h"
#include "RefinementNelderMead.h"
#include "Shouter.h"
#include "Crystal.h"
#include "Diffraction.h"
#include "CSV.h"
#include "Options.h"
#include "Plucker.h"

#define CHECK_DISTANCE_STEP 0.30

Bucket::Bucket()
{
	_mdnode = new MDNode(3);
	_solvBFac = 0;
	_solvScale = 0;
	_wanted = 1;
	
	for (int i = 0; i < 3; i++)
	{
		_averages[i] = 0;
	}
}

void Bucket::scaleSolvent()
{
	if (!_data)
	{
		shout_at_helen("Need diffraction data to scale solvent");
	}
	
	setSolvScale(this, 0.0);
	setSolvBFac(this, 40);
	
	NelderMeadPtr fine;
	fine = NelderMeadPtr(new RefinementNelderMead());
	fine->setJobName("solvent_scale_fine");
	fine->setEvaluationFunction(scaleSolventScore, this);
	fine->setCycles(40);
	fine->addParameter(this, getSolvScale, setSolvScale, 0.4, 0.01, "scale");
	fine->addParameter(this, getSolvBFac, setSolvBFac, 40, 1.0, "bfac");
	fine->setSilent(true);
	fine->refine();

	if (!getCrystal()->isSilent())
	{
		std::cout << "   Solvent B factor: " 
		<< getSolvBFac(this) << std::endl;
	}
	
	/** If we are doing powder analysis we don't actually want
	* 	to add the solvent */
	if (Options::shouldPowder())
	{
		return;
	}	
	
	/** Now add this into the FFT */
	
	FFTPtr fft = getCrystal()->getFFT();
	mat3x3 real2frac = getCrystal()->getReal2Frac();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	FFTPtr fftData = _data->getFFT();
	double nLimit = std::min(fftData->nx, fft->nx);
	nLimit /= 2;

	std::vector<double> fData, fModel;
	CSVPtr csv = CSVPtr(new CSV(2, "fo", "fc"));

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
				double bFacMod = exp(-_solvBFac / four_d_sq);

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
				double bFacMod = exp(-_solvBFac / four_d_sq);

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


void Bucket::applySymOps(CSym::CCP4SPG *spaceGroup)
{
	if (spaceGroup->spg_num == 1)
	{
		return;
	}

	_solvent->applySymmetry(spaceGroup, getCrystal()->isSilent());
}

void Bucket::fourierTransform(int dir)
{
	/* Only care about reciprocal space */
	if (dir == 1)
	{
		CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();

		_solvent->fft(dir);
		applySymOps(spg);
		_solvent->normalise();
		_solvent->data[0][0] = 0;
		_solvent->data[0][1] = 0;
	}
}

Atom *Bucket::nearbyAtom(int index)
{
	if (_atomPtrs.size() <= index)
	{
		return NULL;
	}

	return _atomPtrs[index];
}

bool Bucket::isSolvent(int index)
{
	return (_solvent->data[index][0] > 0.8);
}

bool Bucket::isSolvent(vec3 pos)
{
	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	mat3x3_mult_vec(real2Frac, &pos);
	
	long index = _maskedRegions->elementFromFrac(pos.x, pos.y, pos.z);
	int mask = _maskedRegions->getMask(index);
	
	return (mask == 1);
}

void Bucket::processMaskedRegions()
{
	_maskedRegions = FFTPtr(new FFT(*_solvent));
	_maskedRegions->setupMask();
	FFTPtr fft = getCrystal()->getFFT();
	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	mat3x3 frac2Real = mat3x3_inverse(real2Frac);
	int additions[] = {0, 0, 0};
	double sums[] = {0., 0., 0.};

	std::cout << "Total voxels: " << _maskedRegions->nn << std::endl;

	for (long i = 0; i < _maskedRegions->nn; i++)
	{
		if (_maskedRegions->data[i][0] > 0.8)
		{
			/* solvent */
			_maskedRegions->setMask(i, 1);
			additions[1]++;
			sums[1] += fft->data[i][0];
			continue;
		}
		
		additions[0]++;
		_maskedRegions->setMask(i, 0);
		sums[0] += fft->data[i][0];
		continue;

		vec3 newfrac = _maskedRegions->fracFromElement(i);
		vec3 pos = mat3x3_mult_vec(frac2Real, newfrac);
		
		AtomPtr atom = getCrystal()->getClosestAtom(pos);
		
		if (atom->isHeteroAtom())
		{
			_maskedRegions->setMask(i, 2);
			additions[2]++;
			sums[2] += fft->data[i][0];
			continue;
		}

		additions[0]++;
		sums[0] += fft->data[i][0];
	}
	
	double percentages[3] = {0., 0., 0.};
	
	for (int i = 0; i < 2; i++)
	{
		sums[i] /= (double)additions[i];
		percentages[i] = (double)additions[i] / (double)_maskedRegions->nn;
		_averages[i] = sums[i];
	}
	
	std::cout << std::setprecision(2);
	std::cout << "Protein voxels: " << additions[0];
	std::cout << " (" << percentages[0] * 100 << "%)";
	std::cout << " average value: " << std::setprecision(4)
	<< sums[0] << std::endl;
/*
 * std::cout << "Interface voxels: " << additions[2];
	std::cout << " (" << percentages[1] * 100 << "%)";
	std::cout << " average value: " << std::setprecision(4)
	<< sums[1] << std::endl;
	*/
	std::cout << "Solvent voxels: " << additions[1];
	std::cout << " (" << percentages[1] * 100 << "%)";
	std::cout << " average value: " << std::setprecision(4)
	<< sums[1] << std::endl;
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
void Bucket::populateHistogram(MDNode *node, vec3 centre, vec3 left)
{
	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	vec3 ldiff = vec3_subtract_vec3(left, centre);
	double llength = vec3_length(ldiff);
	double step = CHECK_DISTANCE_STEP;

	FFTPtr fft = getCrystal()->getFFT();

	vec3 transCentre = mat3x3_mult_vec(real2Frac, centre);
	vec3 transLeft = mat3x3_mult_vec(real2Frac, left);
	double centreDensity = fft->getRealFromFrac(transCentre);
	centreDensity -= _averages[_wanted];
	double leftDensity = fft->getRealFromFrac(transLeft);
	leftDensity -= _averages[_wanted];
	
	if (leftDensity < 0) return;
	
	for (double x = -MAX_CHECK_DISTANCE;
	     x <= MAX_CHECK_DISTANCE; x += step)
	{
		for (double y = -MAX_CHECK_DISTANCE;
		     y <= MAX_CHECK_DISTANCE; y += step)
		{
			for (double z = -MAX_CHECK_DISTANCE;
			     z <= MAX_CHECK_DISTANCE; z += step)
			{
				/* difference from centre to right vector */
				vec3 offset = make_vec3(x, y, z);
				double rlength = vec3_length(offset);

				if (rlength < MIN_CHECK_DISTANCE ||
				    rlength > MAX_CHECK_DISTANCE)
				{
					continue;
				}

				vec3 right = vec3_add_vec3(centre, offset);
				
				/* Populate array with result */
				double angle = vec3_angle_with_vec3(ldiff, offset);
				double degrees = rad2deg(angle);
				
				if (degrees < 0) degrees = - degrees;
				
				vec3 transRight = mat3x3_mult_vec(real2Frac, right);
				double rightDensity = fft->getRealFromFrac(transRight);
				rightDensity -= _averages[_wanted];
				
				if (rightDensity < 0) continue;

				double mult = std::min(leftDensity, rightDensity);
				mult = std::min(centreDensity, mult);
				
				vec3 diff = vec3_subtract_vec3(left, offset);
				double dlength = vec3_length(diff);
				
				double dr_ang = vec3_angle_with_vec3(diff, offset);
				double dl_ang = vec3_angle_with_vec3(diff, left);
				
				double dr_deg = rad2deg(dr_ang);
				double dl_deg = rad2deg(dl_ang);
				
				if (mult != mult || degrees != degrees || rlength != rlength)
				{
					continue;	
				}
				
				node->addToNode(mult, 3, llength, rlength, degrees);
				node->addToNode(mult, 3, llength, dlength, dl_deg);
				node->addToNode(mult, 3, rlength, dlength, dr_deg);
			}
		}
	}

}

int Bucket::addAnalysisForSolventPos(MDNode *node, vec3 centre, double distance)
{
	double step = CHECK_DISTANCE_STEP;
	double minDist = distance - step / 3;
	double maxDist = distance + step / 3;

	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	vec3 transCentre = mat3x3_mult_vec(real2Frac, centre);

	double centreDensity = getCrystal()->getFFT()->getRealFromFrac(transCentre);
	centreDensity -= _averages[_wanted];
	
	if (centreDensity < 0) return 0;
	
	std::cout << "Adding data from " << vec3_desc(centre) << std::endl;
	for (double x = -maxDist - step; x <= maxDist + step; x += step)
	{
		for (double y = -maxDist - step; y <= maxDist + step; y += step)
		{
			for (double z = -maxDist - step; z <= maxDist + step; z += step)
			{
				vec3 offset = make_vec3(x, y, z);
				double llength = vec3_length(offset);
				vec3 left = vec3_add_vec3(centre, offset);
				
				/* Temporary removal of some distances */
	
				vec3 ldiff = vec3_subtract_vec3(left, centre);
				if (llength < MIN_CHECK_DISTANCE
				    || llength > MAX_CHECK_DISTANCE)
				{
					continue;
				}

				populateHistogram(node, centre, left);
			}
		}
	}

	CSVPtr csv = CSV::nodeToCSV(node);
	csv->writeToFile("solvent_density_analysis.csv");
	
	return 1;
}

void Bucket::setupNodes(int split)
{
	_mdnode->setDimension(0, MIN_CHECK_DISTANCE, MAX_CHECK_DISTANCE);
	_mdnode->setDimension(1, MIN_CHECK_DISTANCE, MAX_CHECK_DISTANCE);
	_mdnode->setDimension(2, 0, 180);
	_mdnode->splitNode(split, split);
	
}

void Bucket::analyseSolvent(double distance)
{
	const int split = 5;
	setupNodes(split);

	mat3x3 real2Frac = getCrystal()->getReal2Frac();
	mat3x3 frac2Real = mat3x3_inverse(real2Frac);
	
	std::vector<long> shuffled;

	for (long i = 0; i < _maskedRegions->nn; i += 1)
	{
		int val = _maskedRegions->getMask(i);

		if (val != _wanted)
		{
			continue;
		}

		shuffled.push_back(i);
	}
	
	std::random_shuffle(shuffled.begin(), shuffled.end());
	int count = 0;
	
	for (size_t i = 0; i < shuffled.size() && count < 1; i++)
	{
		vec3 newfrac = getCrystal()->getFFT()->fracFromElement(shuffled[i]);
		vec3 centre = mat3x3_mult_vec(frac2Real, newfrac);
		count += addAnalysisForSolventPos(_mdnode, centre, distance);
	}
}

int Bucket::getReallyRandomValues(double left, double *right, double *angle)
{
	*right = (rand() / (double)RAND_MAX) * 6.;
	*angle = (rand() / (double)RAND_MAX) * 180.;
	
	return 1;
}

int Bucket::getRandomValues(double left, double *right, double *angle)
{
	double last = -1;
	Plucker *pluck = NULL;

	for (PluckerItr it = _pluckerMap.begin(); it != _pluckerMap.end(); it++)
	{
		if (last > 0)
		{
			if (it->first < left && it->first > last)
			{
				pluck = it->second;
			}
		}
		
		last = it->first;
	}
	
	if (pluck == NULL)
	{
		return 0;
	}
	
	MDNode *node = static_cast<MDNode *>(pluck->pluck());
	
	*right = node->aveDimension(1);
	*angle = node->aveDimension(2);
	
	return 1;
}
