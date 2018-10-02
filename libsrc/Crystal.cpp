//
//  Crystal.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017-8 Helen Ginn. All rights reserved.
//

#include "Crystal.h"
#include "Bond.h"
#include "fftw3d.h"
#include "vec3.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include "Shouter.h"
#include "Diffraction.h"
#include "Polymer.h"
#include "CSV.h"
#include "FileReader.h"
#include "PDBReader.h"
#include "Atom.h"
#include "RefinementGridSearch.h"
#include "BucketBulkSolvent.h"
#include "BucketMonteCarlo.h"
#include "Options.h"

#include "../libccp4/cmtzlib.h"
#include "../libccp4/csymlib.h"
#include "../libccp4/ccp4_spg.h"
#include "../libccp4/ccp4_general.h"

void Crystal::summary()
{
	std::cout << "|----------------" << std::endl;
	std::cout << "| Crystal summary (" << _filename << "): " << std::endl;
	std::cout << "|----------------" << std::endl;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (i > 0)
		{
			std::cout << "|-------" << std::endl;
		}
		molecule(i)->summary();
	}

	std::cout << "|----------------\n" << std::endl;
}

void Crystal::tieAtomsUp()
{
	if (_tied) return;

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tieAtomsUp();
	}

	_tied = true;
}

void Crystal::addMolecule(MoleculePtr molecule)
{
	if (molecule->getChainID().length() <= 0)
	{
		shout_at_helen("Polymer chain ID is missing while trying\n"\
		               "to interpret PDB file.");
	}

	_molecules[molecule->getChainID()] = molecule;
}

void Crystal::setReal2Frac(mat3x3 mat)
{
	_real2frac = mat;
}

void Crystal::setHKL2Real(mat3x3 mat)
{
	_hkl2real = mat;
}

void Crystal::realSpaceClutter(double maxRes)
{
	if (_fft)
	{
		_fft->setAll(0);	
	}
	
	if (!_fft)
	{
		double sampling = Options::getProteinSampling();

		if (sampling < 0)
		{
			sampling = maxRes / 4.;
			/*
			if (maxRes >= 2.5)
			{
				sampling = maxRes / 5.;
			}
			if (maxRes <= 1.2)
			{
				sampling = maxRes / 3.;
			}
			*/

			Options::setProteinSampling(sampling);
		}
		
		if (sampling >= 1) sampling = 1;

		/* Now create the FFT */
		_fft = FFTPtr(new FFT());
		_difft = FFTPtr(new FFT());

		vec3 uc_dims = empty_vec3();
		vec3 fft_dims = empty_vec3();
		uc_dims.x = mat3x3_length(_hkl2real, 0) / sampling;
		uc_dims.y = mat3x3_length(_hkl2real, 1) / sampling;
		uc_dims.z = mat3x3_length(_hkl2real, 2) / sampling;

		double largest = std::max(uc_dims.x, uc_dims.y);
		largest = std::max(largest, uc_dims.z);

		fft_dims = uc_dims;
		
		if (false)
		{
			fft_dims.x = largest;
			fft_dims.y = largest;
			fft_dims.z = largest;
		}
		
		for (int i = 0; i < 3; i++)
		{
			double val = *(&fft_dims.x + i);
			int fixed = lrint(val);
			if (fixed % 2 > 0) 
			{
				fixed++;
			}
			*(&fft_dims.x + i) = fixed;
		}

		_fft->create(fft_dims.x, fft_dims.y, fft_dims.z);
		_fft->setupMask();

		_difft->create(fft_dims.x, fft_dims.y, fft_dims.z);

		mat3x3 per_voxel = _hkl2real;
		mat3x3_scale(&per_voxel, 1 / fft_dims.x, 
		             1 / fft_dims.y, 1 / fft_dims.z);
		
		_fft->setBasis(per_voxel, 1);
		_difft->setBasis(per_voxel, 1);
	}
	else
	{
		_fft->setAll(0);
		_difft->setAll(0);
	}

	_fft->createFFTWplan(8);
	_difft->createFFTWplan(8);

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->propagateChange();
		molecule(i)->addToMap(_fft, _real2frac);
	}

	if (Options::getAddSolvent())
	{
		if (!_bucket)
		{
			_bucket = BucketPtr(new BucketBulkSolvent());
		}
		
		_bucket->setCrystal(shared_from_this());
		_bucket->addSolvent();
	}
}

double Crystal::totalToScale()
{
	int sum = 0;
	int weighted = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		int weights = 0;
		sum += molecule(i)->totalElectrons(&weights);
		weighted += weights;
	}

	return (sqrt((double)sum / (double)weighted)) * 1.0;
}

void Crystal::omitScan()
{
	int window = 20;
	int half = window / 2;
	
	DiffractionPtr data = Options::getRuntimeOptions()->getActiveData();
	
	for (int j = 0; j < moleculeCount(); j++)
	{
		if (!molecule(j)->isPolymer())
		{
			continue;
		}
		
		PolymerPtr pol = ToPolymerPtr(molecule(j));

		int start = pol->monomerBegin();
		int end = pol->monomerEnd();

		std::cout << "Omit scan starting." << std::endl;

		for (int i = start; i < end; i += half)
		{
			pol->downWeightResidues(i, i + window, 0);
			concludeRefinement(_cycleNum + 1, data);
			pol->downWeightResidues(i, i + window, 1);
		}

		std::cout << "Omit scan complete." << std::endl;
	}
}

void Crystal::reflex()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		PolymerPtr pol = ToPolymerPtr(molecule(i));
		pol->reflex();
	}
}

void Crystal::writeMillersToFile(DiffractionPtr data, std::string prefix)
{
	std::vector<double> bins, ampAves;
	generateResolutionBins(0, _maxResolution, 20, &bins);
	ampAves.resize(bins.size());
	
	for (int i = 0; i < bins.size() - 1; i++)
	{
		ampAves[i] = valueWithDiffraction(data, two_dataset_mean,
		                                 false, bins[i], bins[i + 1]);
	}

	std::string outputFileOnly = prefix + "_" + _filename + "_vbond.mtz";
	getFFT()->writeReciprocalToFile(outputFileOnly, _maxResolution, _spaceGroup,
	                                _unitCell, _real2frac, data->getFFT(),
	                                bins, ampAves);
	std::string outputFile = FileReader::addOutputDirectory(outputFileOnly);
	
	_lastMtz = outputFile;
	
	if (_bucket)
	{
		_bucket->writeMillersToFile(prefix, _maxResolution);	
	}
}

double getNLimit(FFTPtr fftData, FFTPtr fftModel, int axis = 0)
{
	double nLimit = std::min(*(&fftData->nx + axis), 
	                         *(&fftModel->nx + axis));

	nLimit = nLimit - ((int)nLimit % 2);
	nLimit /= 2;

	return nLimit;	
}

vec3 getNLimits(FFTPtr data, FFTPtr fftModel)
{
	vec3 lims;
	lims.x = getNLimit(data, fftModel, 0);
	lims.y = getNLimit(data, fftModel, 1);
	lims.z = getNLimit(data, fftModel, 2);
	return lims;
}

typedef struct
{
	double into;
	double intc;
} IntPair;

void Crystal::scaleAndBFactor(DiffractionPtr data, double *scale, 
                              double *bFactor, FFTPtr model)
{
	if (!model)
	{
		model = _fft;
	}

	std::vector<double> bins;
	generateResolutionBins(6, _maxResolution, 20, &bins);
	std::map<int, std::vector<IntPair> > binRatios;

	FFTPtr fftData = data->getFFT();	
	vec3 nLimits = getNLimits(fftData, model);

	for (int i = -nLimits.z; i < nLimits.z; i++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int k = -nLimits.x; k < nLimits.z; k++)
			{
				int _i = 0; int _j = 0; int _k = 0;
				vec3 ijk = make_vec3(i, j, k);
				CSym::ccp4spg_put_in_asu(_spaceGroup, i, j, k, &_i, &_j, &_k);

				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				double data = fftData->getIntensity(_i, _j, _k);
				double calc = model->getIntensity(i, j, k);

				if (data != data || calc != calc) continue;
				
				IntPair pair;
				pair.into = data;
				pair.intc = calc;
				
				for (int b = 0; b < bins.size() - 1; b++)
				{
					if (length > 1 / bins[b] && length < 1 / bins[b + 1])
					{
						binRatios[b].push_back(pair);
					}
				}
			}
		}
	}
				
	CSVPtr csv = CSVPtr(new CSV(2, "data", "model"));
	std::vector<double> xs, ys;

	for (size_t i = 0; i < bins.size() - 1; i++)
	{
		if (binRatios[i].size() == 0)
		{
			continue;
		}

		double nom = 0;
		double den = 0;
		for (size_t j = 0; j < binRatios[i].size(); j++)
		{
			nom += binRatios[i][j].into;
			den += binRatios[i][j].intc;
		}
		
		double ratio = nom / den;
		double length = (1/bins[i] + 1/bins[i + 1]) / 2;

		double res = 1 / length;
		double four_dsq = 4 * res * res;
		double right_exp = 1 / four_dsq;
		double logratio = log(ratio);

//		std::cout << "Pair: " << right_exp << " " << logratio << std::endl;
		csv->addEntry(2, right_exp, logratio);

		xs.push_back(right_exp);
		ys.push_back(logratio);
	}
	
	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = "bfactor_fit_" + i_to_str(_cycleNum);
	plotMap["height"] = "700";
	plotMap["width"] = "1200";
	plotMap["xHeader0"] = "data";
	plotMap["yHeader0"] = "model";

	plotMap["colour0"] = "black";
	plotMap["xTitle0"] = "1 / (4dd)";
	plotMap["yTitle0"] = "ratio";
	plotMap["style0"] = "line";

	csv->setSubDirectory("correlation_plots");
	csv->plotPNG(plotMap);
	double intercept, gradient;
	regression_line(xs, ys, &intercept, &gradient);

	double k = exp(intercept);
	double b = -gradient;
	
	if (b != b) b = 0;
	
	*scale = k;
	*bFactor = b;
}

double Crystal::valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
                                     bool verbose, double lowRes, double highRes)
{
	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes == 0 ? 1 / _maxResolution : 1 / highRes);

	CSVPtr csv = CSVPtr(new CSV(2, "fo" , "fc"));

	FFTPtr fftData = data->getFFT();	
	vec3 nLimits = getNLimits(fftData, _fft);

	for (int i = -nLimits.z; i < nLimits.z; i++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int k = -nLimits.x; k < nLimits.z; k++)
			{
				int _i = 0; int _j = 0; int _k = 0;
				vec3 ijk = make_vec3(i, j, k);
				CSym::ccp4spg_put_in_asu(_spaceGroup, i, j, k, &_i, &_j, &_k);

				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				if (length < minRes || length > maxRes)
				{
					continue;
				}

				double amp1 = sqrt(fftData->getIntensity(_i, _j, _k));
				double amp2 = sqrt(_fft->getIntensity(i, j, k));

				int isFree = (fftData->getMask(_i, _j, _k) == 0);

				if (amp1 != amp1 || amp2 != amp2)
				{
					continue;
				}

				csv->addEntry(2, amp1, amp2);

				if (!isFree)
				{
					set1.push_back(amp1);
					set2.push_back(amp2);
				}
				else
				{
					free1.push_back(amp1);
					free2.push_back(amp2);
				}
			}
		}
	}

	if (op == r_factor)
	{
		_correlPlotNum++;
		std::string correlName = "correlplot_" + i_to_str(_correlPlotNum);
		csv->setSubDirectory("correlation_plots");
		csv->writeToFile(correlName + ".csv");

		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = correlName;
		plotMap["xHeader0"] = "fo";
		plotMap["yHeader0"] = "fc";
		plotMap["colour0"] = "black";

		plotMap["xTitle0"] = "Fo amplitude";
		plotMap["yTitle0"] = "Fc amplitude";
		plotMap["style0"] = "scatter";
		csv->plotPNG(plotMap);
	}

	_rWork = (*op)(set1, set2);

	if (verbose)
	{
		_ccWork = correlation(set1, set2);
		_ccFree = correlation(free1, free2);
		_rFree = (*op)(free1, free2);
		double diff = _rFree - _rWork; 

		std::cout << "CCwork/CCfree: " << _ccWork * 100 << ", " << _ccFree * 100
		<< " %." << std::endl;

		std::cout << "Rwork/Rfree: " << std::setprecision(4)
		<< _rWork * 100;
		std::cout << ", " << _rFree * 100 <<
		" % (diff: " << diff * 100 << " %)"<<  std::endl;
	}

	return _rWork;
}

void Crystal::applyScaleFactor(double scale, double lowRes, double highRes,
                               double bFactor)
{
	double xLimit = _fft->nx / 2;
	double yLimit = _fft->ny / 2;
	double zLimit = _fft->nz / 2;

	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes <= 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes <= 0 ? FLT_MAX : 1 / highRes);

	/* symmetry issues */
	for (int i = -xLimit; i < xLimit; i++)
	{
		for (int j = -yLimit; j < yLimit; j++)
		{
			for (int k = -zLimit; k < zLimit; k++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);
				long element = _fft->element(i, j, k);

				if (length < minRes || length > maxRes)
				{
					continue;
				}

				double real = _fft->getReal(element);
				double imag = _fft->getImaginary(element);

				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(- bFactor / four_d_sq);
				
				if (i == 0 && j == 0 && k == 0)
				{
					bFacMod = 1;
				}
				
				real *= scale * bFacMod;
				imag *= scale * bFacMod;
				
				if (!std::isfinite(real * real + imag * imag))
				{
					std::cout << "Warning for " << i << " " << j
					<< " " << k << std::endl;
					std::cout << "Scale: " << scale << std::endl;
					std::cout << "bFactor: " << bFacMod << std::endl;
					std::cout << "d: " << d << std::endl;
				}

				_fft->setElement(element, real, imag);
			}
		}
	}
}

void Crystal::scaleSolvent(DiffractionPtr data)
{
	if (!Options::getAddSolvent() || !_bucket)
	{
		return;
	}
	
	_bucket->setData(data);
	_bucket->scaleSolvent();
}

void Crystal::multiplyMap(double scale)
{
	double current = _fft->averageAll();
	std::cout << "Current average: " << current << std::endl;
	_fft->multiplyAll(scale);
	_difft->multiplyAll(scale);
	current = _fft->averageAll();
	std::cout << "New average: " << current << std::endl;
}

double Crystal::getMaxResolution(DiffractionPtr data)
{
	if (_maxResolution <= 0)
	{
		_maxResolution = Options::maxRes();
	}

	if (_maxResolution <= 0)
	{
		if (data)
		{
			_maxResolution = data->getMaxResolution();
		}
		else
		{
			_maxResolution = 1.8;
		}
		std::cout << std::setprecision(2);
	}
	
	return _maxResolution;
}

/* bigger number is more detail */
double Crystal::getMaximumDStar(DiffractionPtr data)
{
	double maxRes = getMaxResolution(data);
	maxRes = 1 / maxRes;
	maxRes *= 2.1;

	return maxRes;
}

double Crystal::getAdjustBFactor()
{
	std::cout << "Adjusting B factor..." << std::endl;
	double change = _bFacFit;
	_bFacFit = 0;
	_realBFactor += change;
	return change;
}

void Crystal::scaleToDiffraction(DiffractionPtr data, bool full)
{
	getMaxResolution(data);
	
	/* First, apply a scale factor to the entire range */
	double totalFc = totalToScale();
	double ratio = valueWithDiffraction(data, &scale_factor_by_sum, false,
	                                    0, _maxResolution);
	applyScaleFactor(totalFc / ratio, 0, 0);

	if (!full)
	{
		/* If non-full scaling has been requested, just an absolute
		 * 	scaling was all that was required. */
		return;
	}
	
	/* If full scaling requested, take global default. */
	ScalingType scaleType = Options::getScalingType();
	double scale, bFactor;
	scaleAndBFactor(data, &scale, &bFactor);
	_bFacFit = bFactor;

	if (scaleType == ScalingTypeAbs)
	{
		/* Same as above, nothing left to do */
		return;
	}
	else if (scaleType == ScalingTypeAbsBFactor)
	{
		std::cout << "Absolute scale: " << scale << " and global B factor: ";
		std::cout << bFactor << std::endl;
		
		applyScaleFactor(totalFc * scale, 0, 0, bFactor);
		
	}
	else if (scaleType == ScalingTypeShell)
	{
		/* Then apply to individual resolution bins */
		std::vector<double> bins;
		generateResolutionBins(0, _maxResolution, 20, &bins);

		/* Extend the final bin by a little bit, so as not to lose any
		 * stragglers. */
		bins[bins.size() - 1] *= 0.95;

		for (int i = 0; i < bins.size() - 1; i++)
		{
			double ratio = valueWithDiffraction(data, &scale_factor_by_sum, false,
			                                    bins[i], bins[i + 1]);
			double scale = totalFc / ratio;

			if (scale != scale)
			{
				scale = 0;
			}

			applyScaleFactor(scale, bins[i], bins[i + 1]);
		}
	}
	else
	{
		std::cout << "Unimplemented scaling method? " << std::endl;
	}
}

void Crystal::scaleComponents(DiffractionPtr data)
{
	scaleToDiffraction(data, false);
	scaleSolvent(data);
	scaleToDiffraction(data);
}

double Crystal::rFactorWithDiffraction(DiffractionPtr data, bool verbose)
{
	double highRes = _maxResolution;
	double lowRes = Options::minRes();

	if (verbose)
	{
		std::cout << "*******************************" << std::endl;
	}

	double rFactor = valueWithDiffraction(data, &r_factor, verbose, 
	                                      lowRes, highRes);

	if (verbose)
	{
		std::cout << "*******************************" << std::endl;
	}

	return rFactor;
}

double Crystal::getDataInformation(DiffractionPtr data, double partsFo,
                                   double partsFc, std::string prefix)
{
	realSpaceClutter(data->getMaxResolution());
	
	std::vector<double> real_calcs;
	int skip = 10;
	for (int i = 0; i < _fft->nn; i += skip)
	{
		real_calcs.push_back(_fft->data[i][0]);
	}
	
	int bad = 0;
	for (int i = 0; i < _fft->nn; i++)
	{
		double val = _fft->data[i][0];
		
		if (val != val)
		{
			bad++;
		}
	}
	
	if (bad > 0)
	{
		std::cout << "There were " << bad << " bad voxels";
		std::cout << " out of " << _fft->nn << "!" << std::endl;
	}
	
	fourierTransform(1, data->getMaxResolution());
	scaleComponents(data);
	
	writeMillersToFile(data, prefix);

	double rFac = rFactorWithDiffraction(data, true);

	FFTPtr fftData = data->getFFT();
	double nLimit = getNLimit(fftData, _fft);
	std::vector<double> set1, set2;

	double lowRes = Options::minRes();
	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = (1 / _maxResolution);
	bool ignoreRfree = Options::ignoreRFree();
	
	if (ignoreRfree)
	{
		warn_user("You should not be ignoring Rfree.\n"\
		          "All your results and any future \n"\
		          "results derived from this are INVALID for\n"\
		          "structure determination.");
	}

	_calcCopy = FFTPtr(new FFT(*_fft));

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = -nLimit; k < nLimit; k++)
			{
				int _h, _k, _l;
				CSym::ccp4spg_put_in_asu(_spaceGroup, i, j, k, &_h, &_k, &_l);

				double amp = sqrt(fftData->getIntensity(_h, _k, _l));
				long index = _fft->element(i, j, k);

				int isAbs = CSym::ccp4spg_is_sysabs(_spaceGroup, i, j, k);
				vec3 ijk = make_vec3(i, j, k);    
				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				bool isRfree = (fftData->getMask(_h, _k, _l) == 0);
				
				if (ignoreRfree)
				{
					isRfree = 0;
				}
				
				if (length < minRes || length > maxRes
				    || isAbs)
				{	
					_fft->setElement(index, 0, 0);
					_difft->setElement(index, 0, 0);

					continue;
				}

				if (amp != amp || isRfree)
				{
					continue;
				}


				vec2 complex;
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImaginary(index);
				double old_amp = sqrt(complex.x * complex.x +
				                      complex.y * complex.y);

				double new_amp = old_amp;
				new_amp = partsFo * amp - partsFc * old_amp;
				new_amp /= old_amp;

				double diff_scale = amp - old_amp;
				diff_scale /= old_amp;

				vec2 diff_complex = complex;

				complex.x *= new_amp;
				complex.y *= new_amp;

				if (amp == amp)
				{
					diff_complex.x *= diff_scale;
					diff_complex.y *= diff_scale;
				}
				
				if (complex.x != complex.x || complex.y != complex.y)
				{
					continue;
				}

				_fft->setElement(index, complex.x, complex.y);
				_difft->setElement(index, diff_complex.x, diff_complex.y);
			}
		}
	}
	
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = -nLimit; k < nLimit; k++)
			{
				bool f000 = (i == 0 && j == 0 && k == 0);
				int _h, _k, _l;
				CSym::ccp4spg_put_in_asu(_spaceGroup, i, j, k, &_h, &_k, &_l);

				double obs_amp = sqrt(fftData->getIntensity(_h, _k, _l));

				vec2 complex;
				int index = _calcCopy->element(i, j, k);

				int isAbs = CSym::ccp4spg_is_sysabs(_spaceGroup, i, j, k);
				vec3 ijk = make_vec3(i, j, k);    
				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				bool isRfree = (fftData->getMask(_h, _k, _l) == 0);

				if (length < minRes || length > maxRes
				    || (isRfree) || isAbs || f000)    
				{	
					continue;
				}

				complex.x = _calcCopy->data[index][0];
				complex.y = _calcCopy->data[index][1];

				double calc_amp = sqrt(_calcCopy->data[index][0] *
				                       _calcCopy->data[index][0] +
				                       _calcCopy->data[index][1] *
				                       _calcCopy->data[index][1]);

				double diff_scale = obs_amp - calc_amp;
				diff_scale /= calc_amp;

				if (diff_scale != diff_scale) continue;

				complex.x *= diff_scale;
				complex.y *= diff_scale;
				
				_difft->data[index][0] = complex.x;
				_difft->data[index][1] = complex.y;
			}
		}
	}
	
	/* Back to real space */
	fourierTransform(-1);
	_difft->fft(-1);

	_calcCopy->scaleToFFT(_fft);

	CSVPtr csv = CSVPtr(new CSV(3, "real_obs", "real_calc", "solvent"));
	std::vector<double> real_mixed, chosen_calc;
	int count = 0;
	
	for (int i = 0; i < _fft->nn; i += skip)
	{
		double obs = _fft->data[i][0];
		double calc = real_calcs[count];
		count++;
		if (calc <= 0) continue;

		double solvent = 0;
		if (_bucket)
		{
			solvent = _bucket->isSolvent(i);
		}

		real_mixed.push_back(obs);
		chosen_calc.push_back(calc);
		csv->addEntry(3, obs, calc, solvent);
	}
	
	double correl = correlation(chosen_calc, real_mixed);
	csv->writeToFile("real_space_cc.csv");
	printf("Real space correlation coefficient: %.3f\n", correl);
	
	if (_bucket)
	{
		_bucket->abandonCalculations();
	}

	return rFac;
}

void Crystal::tiedUpScattering()
{
	double tied = 0;
	double total = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tiedUpScattering(&tied, &total);
		molecule(i)->reportParameters();
	}

	std::cout << std::fixed << std::setprecision(0);
	std::cout << "Tied up " << 100. * sqrt(tied / total) << "% of"\
	" the scattering electrons." << std::endl;
}

// N.B. I say powder, because it reminds me of indexing.
void Crystal::makePowders()
{
	std::cout << "Analysing solvent density." << std::endl;

	if (_bucket)
	{
		_bucket->analyseSolvent(2.0);
	}
}

void Crystal::setAnchors()
{
	std::string anchor = Options::anchorString();
	std::vector<std::string> components = split(anchor, ',');
	
	for (int i = 0; i < components.size(); i++)
	{
		std::string chain = "";
		std::string number = "";

		for (int j = 0; j < components[i].length(); j++)
		{
			if (components[i][j] >= '0' && components[i][j] <= '9')
			{
				number.push_back(components[i][j]);
			}
			else
			{
				chain.push_back(components[i][j]);
			}
		}
		
		int anchorPoint = atoi(number.c_str());
		
		for (int i = 0; i < moleculeCount(); i++)
		{
			std::string thisChain = molecule(i)->getChainID();
			std::string truncated = thisChain.substr(0, chain.length());
			
			if (molecule(i)->isPolymer() && 
			    chain == truncated)
			{
				PolymerPtr polymer = ToPolymerPtr(molecule(i));
				polymer->setAnchor(anchorPoint);
				std::cout << "Setting custom anchor " << anchorPoint;
				std::cout << " for chain " << chain << std::endl;
			}
		}

	}

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() == "Polymer")
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			if (polymer->getAnchor() < 0)
			{
				polymer->findAnchorNearestCentroid();
			}
		}
	}
}

Crystal::Crystal()
{
	_realBFactor = -1;
	_sampleNum = -1;
	_cycleNum = 0;
	_lastRWork = FLT_MAX;
	_bestRWork = FLT_MAX;
	_sinceBestNum = 0;
	_correlPlotNum = 0;
	_tied = false;
	_spaceGroup = NULL;
	_spgNum = 0;
	_spgString = "";
	_maxResolution = 0;
	_solvScale = 0.5;
	_solvBFac = 10;
	_unitCell.resize(6);
	_comments = "";
}

void Crystal::applySymOps(double res)
{
	if (_spaceGroup->spg_num == 1)
	{
		return;
	}

	_fft->applySymmetry(_spaceGroup, res);
}

void Crystal::fourierTransform(int dir, double res)
{
	_fft->fft(dir);

	if (dir == 1)
	{
		applySymOps(res);
	}

	_fft->normalise();
	
	if (_bucket)
	{
		_bucket->fourierTransform(dir, res);
	}
}

void Crystal::vsChangeSampleSize(void *object, double n)
{
	Parser *parser = static_cast<Parser *>(object);
	Crystal *crystal = dynamic_cast<Crystal *>(parser);

	Options::getRuntimeOptions()->setNSamples(n);
}

void Crystal::makePDBs(std::string suffix)
{
	std::vector<std::string> prefices; std::vector<PDBType> pdbTypes;
	prefices.push_back("ensemble_"); pdbTypes.push_back(PDBTypeEnsemble);
	prefices.push_back("average_"); pdbTypes.push_back(PDBTypeAverage);

	std::string path;
	path = FileReader::addOutputDirectory("average_" + suffix + ".pdb");

	_lastAveragePDB = path;

	std::ofstream file;
	file.open(path);
	
	file << PDBReader::writeCryst(_hkl2real, _spaceGroup);

	for (int j = 0; j < moleculeCount(); j++)
	{
		CrystalPtr crystal = shared_from_this();
		file << molecule(j)->makePDB(PDBTypeAverage, crystal); 
	}

	file.close();
	int numConf = 0;
	
	if (moleculeCount() && molecule(0)->atomCount())
	{
		if (!molecule(0)->atom(0)->getModel()->hasExplicitPositions())
		{
			return;
		}

		std::vector<BondSample> *samples; 
		samples = molecule(0)->atom(0)->getExplicitModel()->getManyPositions();
		numConf = samples->size();
	}
	
	path = FileReader::addOutputDirectory("ensemble_" + suffix + ".pdb");
	_lastEnsemblePDB = path;
	std::ofstream ensemble;
	ensemble.open(path);

	for (int j = 0; j < numConf; j++)
	{
		ensemble << "MODEL " << std::setw(8) << j + 1 << std::setw(66)
		<< " " << std::endl;

		for (int i = 0; i < moleculeCount(); i++)
		{
			CrystalPtr crystal = shared_from_this();
			ensemble << molecule(i)->makePDB(PDBTypeEnsemble, crystal, j); 
		}

		ensemble << "TER" << std::setw(80) << " " << std::endl;
		ensemble << "ENDMDL" << std::setw(80) << " " << std::endl;
	}
	
	ensemble.close();
};

void Crystal::writeVagabondFile(int cycleNum)
{
	std::ofstream file;
	std::string filename = "cycle_" + i_to_str(cycleNum) + ".vbond";
	_vbondFile = FileReader::addOutputDirectory(filename);
	file.open(_vbondFile);
	writeToFile(file, 0);
	file.close();
	std::cout << "Written Vagabond model to " << _vbondFile << std::endl;
}

double Crystal::concludeRefinement(int cycleNum, DiffractionPtr data)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->isPolymer())
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			polymer->test();
		}
	}

	_cycleNum = cycleNum;
	std::cout << "*******************************" << std::endl;
	std::cout << "\tCycle " << cycleNum << std::endl;

	std::string refineCount = "cycle_" + i_to_str(cycleNum);
	double rFac = 0;

	if (!data)
	{
		std::cout << "No reflection file has been specified.\n"\
		"Cannot perform map recalculation." << std::endl;
		std::cout << std::setprecision(4);
		realSpaceClutter(1.8);
		_fft->fft(1);
		_fft->normalise();
		_fft->writeReciprocalToFile("calc_" + i_to_str(cycleNum) + ".mtz", 1.8);
		Options::flagDensityChanged();
	}
	else
	{
		rFac = getDataInformation(data, 2, 1, refineCount);
		Options::flagDensityChanged();
	}

	makePDBs(i_to_str(cycleNum));

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() == "Polymer")
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			polymer->graph("bfactor_" + polymer->getChainID() +
			               "_" + i_to_str(cycleNum));
			polymer->closenessSummary();
		}
	}
	
	writeVagabondFile(cycleNum);
	_lastRWork = rFac;
	
	if (rFac < _bestRWork)
	{
		_bestRWork = rFac;
		_sinceBestNum = 0;
	}
	else
	{
		_sinceBestNum++;
	}
	
	differenceAttribution();

	return rFac;
}

void Crystal::setupSymmetry()
{
	mat3x3 hkl2real = mat3x3_from_unit_cell(_unitCell[0], _unitCell[1],
	                                        _unitCell[2], _unitCell[3],
	_unitCell[4], _unitCell[5]);
	mat3x3 real2hkl = mat3x3_inverse(hkl2real);

	setHKL2Real(hkl2real);
	setReal2Frac(real2hkl);
}

std::string Crystal::agreementSummary()
{
	std::ostringstream ss;
	ss << "Rwork/free: " << _rWork * 100 << ", " << _rFree * 100 << "%; ";
	ss << "CCwork/free: " << _ccWork << ", " << _ccFree;
	return ss.str();
}

void Crystal::hydrogenateContents()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		PolymerPtr polymer = ToPolymerPtr(molecule(i));
		
		polymer->hydrogenateContents();
	}
}

void Crystal::fitWholeMolecules(bool translation, bool rotation)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		ToPolymerPtr(molecule(i))->refineAnchorMovements();
	}
}

void Crystal::addProperties()
{
	addStringProperty("filename", &_filename);
	addDoubleProperty("uc_a", &_unitCell[0]);
	addDoubleProperty("uc_b", &_unitCell[1]);
	addDoubleProperty("uc_c", &_unitCell[2]);
	addDoubleProperty("uc_alpha", &_unitCell[3]);
	addDoubleProperty("uc_beta", &_unitCell[4]);
	addDoubleProperty("uc_gamma", &_unitCell[5]);

	addDoubleProperty("r_work", &_rWork);
	addDoubleProperty("r_free", &_rFree);
	addDoubleProperty("cc_work", &_ccWork);
	addDoubleProperty("cc_free", &_ccFree);
	addDoubleProperty("real_b_factor", &_realBFactor);
//	addStringProperty("comments", &_comments);
	addIntProperty("cycles_since_best", &_sinceBestNum);
	addIntProperty("sample_num", &_sampleNum);

	_spgNum = 0;
	_spgString = "";
	if (_spaceGroup)
	{
		_spgNum = _spaceGroup->spg_num;
		_spgString = _spaceGroup->symbol_xHM;
	}

	addIntProperty("spacegroup", &(_spgNum));
	addStringProperty("spg_symbol", &(_spgString));

	for (int i = 0; i < moleculeCount(); i++)
	{
		addChild("molecule", molecule(i));
	}
	
	exposeFunction("recalculate_map", Crystal::vsConcludeRefinement);
	exposeFunction("change_sample_size", Crystal::vsChangeSampleSize);
	exposeFunction("set_shell_scale", Crystal::vsSetShellScale);
	exposeFunction("restore_state", Crystal::vsRestoreState);
}

/* Positive number for ON, negative number for OFF */
void Crystal::vsSetShellScale(void *object, double val)
{
	bool shell = (val > 0);
	Options::getRuntimeOptions()->setScalingType(ScalingTypeShell);
}

void Crystal::vsRestoreState(void *object, double val)
{
	int value = lrint(val);

	Parser *parser = static_cast<Parser *>(object);
	parser->restoreState(value);
}

void Crystal::addObject(ParserPtr object, std::string category)
{
	if (category == "molecule")
	{
		MoleculePtr molecule = ToMoleculePtr(object);
		addMolecule(molecule);
	}
}

void Crystal::postParseTidy()
{
	if (_spgString == "")
	{
		_spaceGroup = CSym::ccp4spg_load_by_ccp4_num(_spgNum);
	}
	else
	{
		_spaceGroup = CSym::ccp4spg_load_by_spgname(_spgString.c_str());
	}

	setupSymmetry();
	_tied = true;
}

AtomPtr Crystal::getClosestAtom(vec3 pos)
{
	AtomPtr atom;
	double small_dist = FLT_MAX;

	CrystalPtr me = shared_from_this();
	for (size_t i = 0; i < moleculeCount(); i++)
	{
		AtomPtr tmp = molecule(i)->getClosestAtom(me, pos);
		vec3 tmp_pos = tmp->getAsymUnitPosition(me);

		vec3 diff = vec3_subtract_vec3(pos, tmp_pos);
		double dist = vec3_length(diff);
		
		if (tmp && dist < small_dist)
		{
			small_dist = dist;
			atom = tmp;
		}
	}
	
	return atom;
}

std::vector<AtomPtr> Crystal::getCloseAtoms(std::vector<AtomPtr> atoms, 
                                            double tol)
{
	std::vector<AtomPtr> extra;
	
	for (size_t i = 0; i < atoms.size(); i++)
	{
		std::vector<AtomPtr> clAtoms = getCloseAtoms(atoms[i], tol);
		
		for (int j = 0; j < clAtoms.size(); j++)
		{
			if (std::find(atoms.begin(), atoms.end(), clAtoms[j]) !=
			    atoms.end())
			{
				continue;
			}

			if (std::find(extra.begin(), extra.end(), clAtoms[j]) !=
			    extra.end())
			{
				continue;
			}

			extra.push_back(clAtoms[j]);
		}
	}
	
	return extra;
}

std::vector<AtomPtr> Crystal::getCloseAtoms(AtomPtr one, double tol, bool cache)
{
	std::vector<AtomPtr> atoms;

	for (int i = 0; i < moleculeCount(); i++)
	{
		std::vector<AtomPtr> someAtoms = molecule(i)->getCloseAtoms(one, tol, cache);
		atoms.reserve(atoms.size() + someAtoms.size());
		atoms.insert(atoms.end(), someAtoms.begin(), someAtoms.end());
	}

	return atoms;
}

vec3 Crystal::snapToGrid(vec3 pos)
{
	mat3x3_mult_vec(_real2frac, &pos);
	pos.x *= _fft->nx;
	pos.y *= _fft->ny;
	pos.z *= _fft->nz;
	
	pos.x = lrint(pos.x);
	pos.y = lrint(pos.y);
	pos.z = lrint(pos.z);
	
	pos.x /= _fft->nx;
	pos.y /= _fft->ny;
	pos.z /= _fft->nz;
	mat3x3_mult_vec(_hkl2real, &pos);
	return pos;
}

std::vector<AtomPtr> Crystal::getAtomsInBox(vec3 target, double tolx,
                                            double toly, double tolz)
{
	std::vector<AtomPtr> atoms;

	for (int i = 0; i < moleculeCount(); i++)
	{
		for (int j = 0; j < molecule(i)->atomCount(); j++)
		{
			AtomPtr anAtom = molecule(i)->atom(j);
			vec3 pos = anAtom->getAbsolutePosition();
			
			vec3 diff = vec3_subtract_vec3(pos, target);
			
			if (fabs(diff.x) > tolx || fabs(diff.y) > toly
			    || fabs(diff.z) > tolz)
			{
				continue;
			}
			
			atoms.push_back(anAtom);
		}
	}

	return atoms;
}

void Crystal::clearCloseCache()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->clearCloseCache();
	}	
}

double Crystal::vsConcludeRefinement(void *object)
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();
	
	Parser *parser = static_cast<Parser *>(object);
	Crystal *crystal = dynamic_cast<Crystal *>(parser);
	crystal->_cycleNum++;
	int num = crystal->_cycleNum;
	crystal->concludeRefinement(num, data);
	crystal->saveState();
	
	Options::getRuntimeOptions()->agreementSummary();
	std::string agreement = crystal->agreementSummary();
	crystal->addComment("Recalculated: " + agreement);
	
	return 0;
}

void Crystal::postRestoreState()
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();
	CrystalPtr crystal = options->getActiveCrystal();
	
	_cycleNum++;
	concludeRefinement(_cycleNum, data);
	crystal->saveState();
}

void Crystal::openInCoot()
{
	std::string command = "coot " + _lastMtz + " " + _lastEnsemblePDB
	+ " " + _lastAveragePDB + " &\n";
	
	std::cout << "Terminal command: " << command << std::endl;

	system(command.c_str());
}

std::vector<AtomPtr> Crystal::getHydrogenBonders()
{
	std::vector<AtomPtr> returns;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->isWaterNetwork()) continue;

		std::vector<AtomPtr> bonders = molecule(i)->getHydrogenBonders();
		
		returns.reserve(returns.size() + bonders.size());
		returns.insert(returns.end(), bonders.begin(), bonders.end());
	}
	
	return returns;
}

int Crystal::getSampleNum()
{
	if (Options::getNSamples() >= 0)
	{
		_sampleNum = Options::getNSamples();
		Options::setNSamples(-1);
	}

	if (_sampleNum < 0) 
	{
		_sampleNum = 120;
	}
	
	return _sampleNum;
}

void Crystal::differenceAttribution()
{
	if (!_bucket)
	{
		return;
	}

	double sum_solvent = 0;
	double num_solvent = 0;
	double sum_all = 0;
	double sum_side = 0;
	double sum_back = 0;
	double sum_hetatm = 0;
	double num_all = 0;
	double num_side = 0;
	double num_back = 0;
	double num_hetatm = 0;
	
	for (int i = 0; i < _difft->nn; i++)
	{
		vec3 frac = _difft->fracFromElement(i);
		
		if (frac.x > _spaceGroup->mapasu_zero[0] ||
		    frac.y > _spaceGroup->mapasu_zero[1] ||
		    frac.z > _spaceGroup->mapasu_zero[2]) 
		{
			continue;
		}
		
		Atom *atom = _bucket->nearbyAtom(i);
		double density = _difft->data[i][0];
		density *= density;
		
		sum_all += fabs(density);
		num_all++;
		
		if (!atom)
		{
			sum_solvent += fabs(density);
			num_solvent++;
		}
		else
		{
			if (atom->isHeteroAtom())
			{
				sum_hetatm += fabs(density);
				num_hetatm++;
			}
			else if (atom->isBackbone() || atom->isBackboneAndSidechain())
			{
				sum_back += fabs(density);
				num_back++;
			}
			else
			{
				sum_side += fabs(density);
				num_side++;
			}
		}
	}
	
	double num_atom = (num_all - num_solvent);
	
	double sum_protein = sum_all - sum_solvent;
	double nearAtom = sum_protein / sum_all * 100;
	double inSolvent = sum_solvent / sum_all * 100;
	
	double cNearAtom = (nearAtom / num_atom);
	double cSolvent = (inSolvent / num_solvent);
	double cAll = cNearAtom + cSolvent;
	
	double backbone = sum_back / sum_protein * 100;
	double sidechain = sum_side / sum_protein * 100;
	double hetatm = sum_hetatm / sum_protein * 100;
	double cBackbone = sum_back / num_back;
	double cSidechain = sum_side / num_side;
	double cHetatm = sum_hetatm / num_hetatm;
	double cProtein = cHetatm + cSidechain + cBackbone;
	
	cNearAtom *= 100 / cAll;
	cSolvent *= 100 / cAll;
	
	cBackbone *= 100 / cProtein;
	cSidechain *= 100 / cProtein;
	cHetatm *= 100 / cProtein;
	
	std::cout << std::setprecision(2);
	std::cout << "---------------------------------------------------" 
	<< std::endl;
	std::cout << "| Difference density locations  | Volume-corrected"
	<< std::endl;
	std::cout << "---------------------------------------------------" 
	<< std::endl;
	std::cout << "| Near modelled atoms: " << std::setw(5) << nearAtom << "%   ";
	std::cout << "|  " << std::setw(5) << cNearAtom << "%" << std::endl;
	std::cout << "| Far from any atom:   " << std::setw(5) << inSolvent << "%   ";
	std::cout << "|  " << std::setw(5) << cSolvent << "%" << std::endl;
	std::cout << "---------------------------------------------------" 
	<< std::endl;
	std::cout << "| Protein breakdown: " << std::endl;
	std::cout << "| Backbone atoms:      " << std::setw(5) << backbone << "%   ";
	std::cout << "|  " << std::setw(5) << cBackbone << "%" << std::endl;
	std::cout << "| Sidechain atoms:     " << std::setw(5) << sidechain << "%   ";
	std::cout << "|  " << std::setw(5) << cSidechain << "%" << std::endl;
	std::cout << "| HETATM atoms:        " << std::setw(5) << hetatm << "%   ";
	std::cout << "|  " << std::setw(5) << cHetatm << "%" << std::endl;
	std::cout << "---------------------------------------------------" 
	<< std::endl;
	std::cout << std::endl;
}

double Crystal::getRealBFactor()
{
	if (Options::getGlobalBFactor() >= 0)
	{
		_realBFactor = Options::getGlobalBFactor();
		Options::resetGlobalBFactor();
	}
	
	if (_realBFactor <= 0)
	{
		_realBFactor = 20;
	}
	
	return _realBFactor;
}

