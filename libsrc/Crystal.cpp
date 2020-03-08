// Vagabond
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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <time.h>

#include "Crystal.h"
#include "FFT.h"
#include "Bond.h"
#include "vec3.h"
#include "Shouter.h"
#include "Diffraction.h"
#include "Polymer.h"
#include "CSV.h"
#include "FileReader.h"
#include "PDBReader.h"
#include "Atom.h"
#include "Anchor.h"
#include "FlexGlobal.h"
#include "RefinementGridSearch.h"
#include "RefinementNelderMead.h"
#include "Bucket.h"
#include "Options.h"
#include "WeightedMap.h"

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
	if (!molecule)
	{
		return;
	}

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

void Crystal::refinePolymers(RefinementType type)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		ToPolymerPtr(molecule(i))->refine(shared_from_this(), type);
	}

}

void Crystal::refineCrude()
{
	refinePolymers(RefinementCrude);
}

void Crystal::refinePositions()
{
	refinePolymers(RefinementModelPos);
}

void Crystal::refineSidechains()
{
	refinePolymers(RefinementSidechain);
}

void Crystal::refineSidechainPositions()
{
	refinePolymers(RefinementSidePos);
}

bool Crystal::refineIntraMovements(bool magic)
{
	bool changed = false;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		changed |= ToPolymerPtr(molecule(i))->refineLocalFlexibility(magic);
	}
	
	return changed;
}

double Crystal::getProteinSampling()
{
	double sampling = Options::getProteinSampling();
	double maxRes = getMaxResolution(_data);
	
	if (sampling < 0)
	{
		sampling = maxRes / 3.;

		if (sampling >= 1.0)
		{
			sampling = 1.0;
		}

		Options::setProteinSampling(sampling);
	}

	return sampling;
}

void Crystal::realSpaceClutter(double maxRes)
{
	if (!_fft)
	{
		double sampling = getProteinSampling();
		
		/* work out what our nx/ny/nz dimensions should be */
		vec3 uc_dims = empty_vec3();
		vec3 fft_dims = empty_vec3();
		uc_dims.x = mat3x3_length(_hkl2real, 0) / sampling;
		uc_dims.y = mat3x3_length(_hkl2real, 1) / sampling;
		uc_dims.z = mat3x3_length(_hkl2real, 2) / sampling;

		fft_dims = uc_dims;
		
		/* even-ise them */
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

		/* Now create the FFTs */
		_fft = VagFFTPtr(new VagFFT(fft_dims.x, fft_dims.y, fft_dims.z));
		_difft = VagFFTPtr(new VagFFT(fft_dims.x, fft_dims.y, fft_dims.z));

//		_fft->setupMask();

		if (!_original)
		{
			_original = VagFFTPtr(new VagFFT(fft_dims.x, 
			                                 fft_dims.y, fft_dims.z));
			prepareFFT(_original);
			setupOriginalMap();
		}

		prepareFFT(_fft);
		prepareFFT(_difft);
	}
	else
	{
		_fft->wipe();
		_difft->wipe();
		_fft->setStatus(FFTRealSpace);
		_difft->setStatus(FFTRealSpace);
	}

	refreshAnchors();
	refreshPositions();

	addToMap(_fft);
}

void Crystal::prepareFFT(VagFFTPtr ft)
{
	ft->setSpaceGroup(_spaceGroup);
	ft->setUnitCell(_unitCell);
	ft->makePlans();
	ft->setStatus(FFTRealSpace);
}

void Crystal::setupOriginalMap()
{
	if (!_data)
	{
		return;
	}
	
	if (!_data->getOriginal())
	{
		return;
	}

	_original = VagFFTPtr(new VagFFT(*_data->getOriginal()));
	_original->setStatus(FFTReciprocalSpace);
	_original->setSpaceGroup(_spaceGroup);
	_original->setUnitCell(_unitCell);
	_original->makePlans();

	/* to real space */
	_original->fft(FFTReciprocalToReal);
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
	
	double unmodelled = Options::getUnmodelledFraction();
	
	double extra = (1 / (1 - unmodelled));
	
	if (unmodelled > 0)
	{
		std::cout << "Accounting for unmodelled percentage: " <<
		unmodelled  * 100 << std::endl;
		std::cout << "Total content: " <<
		extra * 100 << "%" << std::endl;
	}

	return (sqrt((double)sum / (double)weighted)) * extra;
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

void Crystal::writeMillersToFile(DiffractionPtr data, std::string prefix)
{
	if (_silent)
	{
		return;
	}

	std::string outputFileOnly = prefix + "_" + _filename + "_vbond.mtz";
	std::string outputFile = FileReader::addOutputDirectory(outputFileOnly);
	
	_lastMtz = outputFile;
	
	if (_bucket)
	{
		_bucket->writeMillersToFile(prefix, _maxResolution);	
	}
}

typedef struct
{
	double into;
	double intc;
} IntPair;

void Crystal::scaleAndBFactor(DiffractionPtr data, double *scale, 
                              double *bFactor, VagFFTPtr model)
{
	if (!model)
	{
		model = _fft;
	}

	std::vector<double> bins;
	generateResolutionBins(6, _maxResolution, 20, &bins);
	std::map<int, std::vector<IntPair> > binRatios;

	VagFFTPtr fftData = data->getFFT();	
	vec3 nLimits = getNLimits(fftData, model);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				int _i = 0; int _j = 0; int _k = 0;
				vec3 ijk = make_vec3(i, j, k);
				CSym::ccp4spg_put_in_asu(_spaceGroup, i, j, k, &_i, &_j, &_k);

				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				double data = fftData->getReal(_i, _j, _k);
				double calc = model->getAmplitude(i, j, k);

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
				
	CSVPtr csv = CSVPtr(new CSV(3, "res", "data", "model"));
	std::vector<double> xs, ys, zs;
	_resBinAves = std::map<double, double>();

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
		
		nom /= (double)binRatios[i].size();
		den /= (double)binRatios[i].size();
		
		/* Get average F-obs per resolution shell */
		_resBinAves[bins[i]] = sqrt(nom);
		
		double length = (1/bins[i] + 1/bins[i + 1]) / 2;

		double res = 1 / length;
		double four_dsq = 4 * res * res;
		double right_exp = 1 / four_dsq;
		double logcalc = log(den);
		double logobs = log(nom);

		csv->addEntry(3, right_exp, logobs, logcalc);
		xs.push_back(right_exp);
		ys.push_back(logobs);
		zs.push_back(logcalc);
	}
	
	if (true || Options::makeDiagnostics())
	{
		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = "bfactor_fit_" + i_to_str(_cycleNum);
		plotMap["height"] = "700";
		plotMap["width"] = "1200";
		plotMap["xHeader0"] = "res";
		plotMap["xHeader1"] = "res";
		plotMap["yHeader0"] = "data";
		plotMap["yHeader1"] = "model";

		plotMap["colour0"] = "blue";
		plotMap["colour1"] = "black";
		plotMap["xTitle0"] = "1 / (4dd)";
		plotMap["yTitle0"] = "log(val)";
		plotMap["style0"] = "line";
		plotMap["style1"] = "line";

		csv->setSubDirectory("correlation_plots");
		csv->plotPNG(plotMap);
	}

	double intercept, gradient;
	regression_line(xs, ys, &intercept, &gradient);
	if (!_silent)
	{
		std::cout << "Wilson plot  (data): " << -gradient << std::endl;
	}
	regression_line(xs, zs, &intercept, &gradient);
	if (!_silent)
	{
		std::cout << "Wilson plot (model): " << -gradient << std::endl;
	}

	double k = exp(intercept);
	double b = -gradient;
	
	if (b != b) b = 0;
	
	*scale = k;
	*bFactor = b;
}

void Crystal::applyShellFactors(DiffractionPtr data)
{
	VagFFTPtr fftData = data->getFFT();	
	vec3 nLimits = getNLimits(fftData, _fft);
	mat3x3 tmp = mat3x3_transpose(_real2frac);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				int _i = 0; int _j = 0; int _k = 0;
				vec3 ijk = make_vec3(i, j, k);
				long element = _fft->element(i, j, k);

				mat3x3_mult_vec(tmp, &ijk);
				double length = vec3_length(ijk);
				double real_space = 1 / length;

				int index = -1;
				
				for (int l = 0; l < _shells.size(); l++)
				{
					if (real_space <= _shells[l].minRes &&
					    real_space > _shells[l].maxRes)
					{
						index = l;
						break;
					}
				}

				if (index < 0)
				{
					continue;
				}

				double real = _fft->getReal(element);
				double imag = _fft->getImag(element);
				
				double scale = _shells[index].scale;
				
				if (scale != scale || !std::isfinite(scale))
				{
					scale = 0;
				}

				real /= scale;
				imag /= scale;

				_fft->setElement(element, real, imag);
			}
		}
	}

}

double Crystal::valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
                                     bool allShells, bool verbose)
{
	std::vector<double> set1, set2, free1, free2;

	CSVPtr csv = CSVPtr(new CSV(5, "h", "k", "l", "fo" , "fc"));

	VagFFTPtr fftData = data->getFFT();	
	vec3 nLimits = getNLimits(fftData, _fft);
	mat3x3 tmp = mat3x3_transpose(_real2frac);
	
	double maxRes = 0;
	double minRes = Options::minRes();

	if (_shells.size())
	{
		maxRes = _shells.back().maxRes;
	}

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				int _i = 0; int _j = 0; int _k = 0;
				vec3 ijk = make_vec3(i, j, k);
				CSym::ccp4spg_put_in_asu(_spaceGroup, i, j, k, &_i, &_j, &_k);

				mat3x3_mult_vec(tmp, &ijk);
				double length = vec3_length(ijk);
				double real = 1 / length;
				
				if (real < maxRes || real > minRes)
				{
					continue;
				}

				int index = -1;
				
				for (int l = 0; l < _shells.size() && !allShells; l++)
				{
					if (real <= _shells[l].minRes &&
					    real > _shells[l].maxRes)
					{
						index = l;
						break;
					}
				}

				if (index < 0 && !allShells)
				{
					continue;
				}

				long fi = fftData->element(_i, _j, _k);
				double amp1 = fftData->getReal(fi);
				double amp2 = _fft->getAmplitude(i, j, k);

				float rfree = fftData->getScratchComponent(fi, 0, 0);
				int isFree = (rfree < 0.5);

				if (amp1 != amp1 || amp2 != amp2 || rfree != rfree)
				{
					continue;
				}

				if (!allShells)
				{
					if (!isFree)
					{
						_shells[index].work1.push_back(amp1);
						_shells[index].work2.push_back(amp2);
					}
					else
					{
						_shells[index].free1.push_back(amp1);
						_shells[index].free2.push_back(amp2);
					}
				}
				else
				{
					csv->addEntry(3, (double)i, (double)j, (double)k, 
					              amp1, amp2);

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
	}
	
	if (op == r_factor && allShells)
	{
		_correlPlotNum++;
		std::string correlName = "correlplot_" + i_to_str(_correlPlotNum);
		csv->setSubDirectory("correlation_plots");
		
		if (Options::makeDiagnostics())
		{
			csv->writeToFile(correlName + ".csv");
		}

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

	if (allShells && !_silent)
	{
		_rWork = (*op)(set1, set2);
		_ccWork = correlation(set1, set2);
		_ccFree = correlation(free1, free2);
		_rFree = (*op)(free1, free2);
		double diff = _rFree - _rWork; 

		if (verbose)
		{
			std::cout << "CCwork/CCfree: " << _ccWork * 100 << ", " << _ccFree * 100
			<< " %." << std::endl;

			std::cout << "Rwork/Rfree: " << std::setprecision(4)
			<< _rWork * 100;
			std::cout << ", " << _rFree * 100 <<
			" % (diff: " << diff * 100 << " %)"<<  std::endl;
		}
	}
	
	if (!allShells)
	{
		for (int i = 0; i < _shells.size(); i++)
		{
			_shells[i].scale = scale_factor_by_sum(_shells[i].work1,
			                                       _shells[i].work2);
		}
	}

	return _rWork;
}

void Crystal::applyScaleFactor(double scale, double lowRes, double highRes,
                               double bFactor)
{
	double xLimit = _fft->nx() / 2;
	double yLimit = _fft->ny() / 2;
	double zLimit = _fft->nz() / 2;

	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes <= 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes <= 0 ? FLT_MAX : 1 / highRes);

	mat3x3 tmp = mat3x3_transpose(_real2frac);

	/* symmetry issues */
	for (int i = -xLimit; i < xLimit; i++)
	{
		for (int j = -yLimit; j < yLimit; j++)
		{
			for (int k = -zLimit; k < zLimit; k++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(tmp, &ijk);
				double length = vec3_length(ijk);
				long element = _fft->element(i, j, k);

				if (length < minRes || length > maxRes)
				{
					continue;
				}

				double real = _fft->getReal(element);
				double imag = _fft->getImag(element);

				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(- bFactor / four_d_sq);
				
				if (i == 0 && j == 0 && k == 0)
				{
					bFacMod = 1;
				}
				
				real *= scale * bFacMod;
				imag *= scale * bFacMod;

				_fft->setElement(element, real, imag);
			}
		}
	}
}

bool Crystal::undoIfWorse()
{
	if (_lastRWork < _rWork)
	{
		std::cout << "Decided to go back one state "
		" results (Rwork = " << _lastRWork * 100 << "%, now " 
		<< _rWork * 100 << "%) "
		"due to Rwork rise since." << std::endl;

		restoreState(-1);
		_sinceBestNum = 0;
		return true;
	}
	else
	{
		std::cout << "Rwork has improved since last cycle" << std::endl;
		return false;
	}
}

bool Crystal::returnToBestState()
{
	if (_rWork > _bestRWork)
	{
		std::cout << "Decided to undo to state " << _bestState << 
		" results (Rwork = " << _bestRWork * 100 << "%, now " 
		<< _rWork * 100 << "%) "
		"due to Rwork rise since." << std::endl;

		restoreState(_bestState);
		_sinceBestNum = 0;
		_bestState = stateCount() - 1;
		return true;
	}
	else
	{
		std::cout << "Not undoing result due to improved Rwork" << std::endl;
		return false;
	}
}

void Crystal::scaleSolvent(DiffractionPtr data)
{
	if (!Options::getAddSolvent())
	{
		return;
	}

	_bucket = Bucket::chosenBucket();
		
	if (!_bucket)
	{
		return;
	}

	_bucket->setCrystal(shared_from_this());
	_bucket->addSolvent();

	std::cout << "Ave calculated density before: " << _fft->averageAll() << 
	" electrons/A^(-3)" << std::endl;
	std::cout << "Ave partial density before: " << 
	_bucket->getSolvent()->averageAll() << 
	" electrons/A^(-3)" << std::endl;

	_bucket->fourierTransform(1);
	_bucket->setData(data);
	_bucket->scalePartialStructure();
	_bucket->reportScale();
	_bucket->postScaleWork();
}

double Crystal::getMaxResolution(DiffractionPtr data)
{
	if (!data)
	{
		data = _data;
	}

	std::cout << std::setprecision(2);
	if (_maxResolution <= 0)
	{
		_maxResolution = Options::maxRes();
	}

	if (_maxResolution <= 0)
	{
		if (data)
		{
			_maxResolution = data->getMaxResolution();
			std::cout << "Maximum resolution: " << _maxResolution <<
			" Å." << std::endl;
		}
		else
		{
			_maxResolution = 1.8;
		}

	}
	
	return _maxResolution;
}

double Crystal::getAdjustBFactor()
{
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
	double ratio = valueWithDiffraction(data, &scale_factor_by_sum, 
	                                    true, false);

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

	/* We want to know the Wilson plot results */
	scaleAndBFactor(data, &scale, &bFactor);
	_bFacFit = bFactor;

	if (scaleType == ScalingTypeAbs)
	{
		/* Same as above, nothing left to do */
		/* just for vagadensity */
		makeShells(&_shells, 0, _maxResolution);
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
		double minRes = Options::minRes();
		makeShells(&_shells, minRes, _maxResolution);

		valueWithDiffraction(data, &scale_factor_by_sum, false, false);

		applyShellFactors(data);
	}
	else
	{
		std::cout << "Unimplemented scaling method? " << std::endl;
	}
}

void Crystal::scaleAnyPartialSet()
{
	if (!Options::usePartial())
	{
		return;
	}

	PartialStructurePtr structure;
	structure = PartialStructurePtr(new PartialStructure());
	
	VagFFTPtr partial = _data->getPartial();
	structure->setData(_data);
	structure->setCrystal(shared_from_this());
	
	if (!partial)
	{
		std::string warn = "Asked to scale partial data set, but no\n"
		"partial data set found in FPART/PHIPART columns.";
		warn_user(warn);
		return;
	}
	
	structure->setStructure(partial);
	structure->scalePartialStructure();
	structure->reportScale();
}

void Crystal::scaleComponents(DiffractionPtr data)
{
	std::cout << "Scaling model to data..." << std::endl;
	/* Just scale using an absolute value only */
	scaleToDiffraction(data, false);
	scaleSolvent(data);
	scaleAnyPartialSet();
	/* Scale using the favoured mechanism (e.g. per-shell) */
	scaleToDiffraction(data);
}

double Crystal::rFactorWithDiffraction(DiffractionPtr data, bool verbose)
{
	double highRes = _maxResolution;
	double lowRes = Options::minRes();

	if (verbose && !_silent)
	{
		std::cout << "*******************************" << std::endl;
	}

	double rFactor = valueWithDiffraction(data, &r_factor, true, true);


	if (verbose && !_silent)
	{
		std::cout << "*******************************" << std::endl;
	}

	return rFactor;
}

double Crystal::getDataInformation(DiffractionPtr data, double partsFo,
                                   double partsFc, std::string prefix)
{
	realSpaceClutter(data->getMaxResolution());
	
	fourierTransform(1);
	
	scaleComponents(data); /* change name of function to write solvent only */
	
	writeMillersToFile(data, prefix);

	double rFac = rFactorWithDiffraction(data, true);
	
	WeightedMap wMap;
	wMap.setCrystalAndData(shared_from_this(), _data);
	wMap.createWeightedMaps();

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
	int flex = 0;
	int pos = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tiedUpScattering(&tied, &total);
		molecule(i)->reportParameters();
		molecule(i)->addParamCounts(&pos, &flex);
	}
	
	std::cout << "Total positional params: " << pos << std::endl;
	std::cout << "Total flexibility params: " << flex << std::endl;
	std::cout << "Total params: " << pos + flex << std::endl;

	std::cout << std::fixed << std::setprecision(0);
	std::cout << "Tied up " << 100. * sqrt(tied / total) << "% of"\
	" the scattering electrons." << std::endl;
	std::cout << std::endl;
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
	_lastLocalCC = 0;
	_bestState = 0;
	_probeRadius = 0.;
	_silent = false;
	_largestNum = -INT_MAX;
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
	_unitCell.resize(6);
	_comments = "";
}

void Crystal::applySymOps()
{
	if (_spaceGroup->spg_num == 1)
	{
		return;
	}

	std::cout << "Protein: ";
	_fft->applySymmetry(false);
}

void Crystal::fourierTransform(int dir, double res)
{
	if (dir == 1)
	{
		_fft->fft(FFTRealToReciprocal);
		applySymOps();
	}
	else
	{
		_fft->fft(FFTReciprocalToReal);
	}
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
	_data = data;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->isPolymer())
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			polymer->test();
		}
	}

	_cycleNum = cycleNum;
	
	if (!_silent)
	{
		writeVagabondFile(cycleNum);
	}
	
	if (!_silent)
	{
		std::cout << "*******************************" << std::endl;
		std::cout << "\tCycle " << cycleNum << std::endl;
	}

	std::string refineCount = "cycle_" + i_to_str(cycleNum);
	double rFac = 0;

	if (!data)
	{
		std::cout << "No reflection file has been specified.\n"\
		"Cannot perform map recalculation." << std::endl;
		std::cout << std::setprecision(4);
		realSpaceClutter(1.8);
		_fft->fft(FFTRealToReciprocal);
//		_fft->writeReciprocalToFile("calc_" + i_to_str(cycleNum) + ".mtz", 1.8);
		_fft->writeToFile("calc_" + i_to_str(cycleNum) + ".mtz", 1.8);
		Options::flagDensityChanged();
	}
	else
	{
		rFac = getDataInformation(data, 2, 1, refineCount);
		Options::flagDensityChanged();
	}
	
	if (_silent)
	{
		return 0;
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

MotionPtr Crystal::getOverallMotion()
{
	if (_motions.size())
	{
		for (int i = 0; i < motionCount(); i++)
		{
			if (_motions[i]->getName() == "all")
			{
				return _motions[i];
			}
		}
	}

	return MotionPtr();
}

void Crystal::makeOverallMotion()
{
	if (!_tied)
	{
		return;
	}

	if (getOverallMotion() != MotionPtr())
	{
		return;
	}
	
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer() ||
		    !molecule(i)->isTied())
		{
			continue;
		}

		PolymerPtr pol = ToPolymerPtr(molecule(i));
		pol->getAnchorModel()->atLeastOneMotion();
	}
}

size_t Crystal::polymerCount()
{
	size_t count = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		count++;
	}

	return count;
}

void Crystal::fitWholeMolecules(bool recip)
{
	std::cout << "Refining in " << 
	(recip ? "reciprocal " : "real ") << 
	"space." << std::endl;

	for (int i = 0; i < _motions.size(); i++)
	{
		if (_motions[i]->getName() == "all")
		{
			continue;
		}

		_motions[i]->refine(recip);
	}
}

void Crystal::rigidBodyRefinement()
{
	for (int i = 0; i < motionCount(); i++)
	{
		_motions[i]->rigidRefine();
	}

}

bool Crystal::calibrateAllMolecules()
{
	FlexGlobal target;
	AtomGroupPtr all = AtomGroupPtr(new AtomGroup());
	NelderMeadPtr neld;
	neld = NelderMeadPtr(new RefinementNelderMead());
	int count = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}

		count++;
		PolymerPtr pol = ToPolymerPtr(molecule(i));
		AtomGroupPtr allBackbone = pol->getAllBackbone();
		all->addAtomsFrom(allBackbone);
		AnchorPtr anch = pol->getAnchorModel();

		neld->addParameter(&*anch, Anchor::sgetBFactor, Anchor::ssetBFactor, 
						   8.0, 0.01, "b_" + pol->getChainID());
	}

	if (count <= 1)
	{
		return false;
	}
	
	target.setCrystal(shared_from_this());
	target.setAtomGroup(all);

	neld->setEvaluationFunction(FlexGlobal::score, &target);
	neld->setVerbose(true);
	neld->setCycles(60);
	FlexGlobal::score(&target);
	
	all->scoreWithMap(ScoreTypeCorrel, shared_from_this(), "all_before");

	neld->refine();

	all->scoreWithMap(ScoreTypeCorrel, shared_from_this(), "all_after");
	
	return true;
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
	addDoubleProperty("probe_radius", &_probeRadius);
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
	
	for (size_t i = 0; i < motionCount(); i++)
	{
		addChild("motion", _motions[i]);
	}
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

	if (category == "motion")
	{
		MotionPtr motion = ToMotionPtr(object);
		addMotion(motion);
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
		vec3 tmp_pos = tmp->getPositionInAsu();

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
	pos.x *= _fft->nx();
	pos.y *= _fft->ny();
	pos.z *= _fft->nz();
	
	pos.x = lrint(pos.x);
	pos.y = lrint(pos.y);
	pos.z = lrint(pos.z);
	
	pos.x /= _fft->nx();
	pos.y /= _fft->ny();
	pos.z /= _fft->nz();
	mat3x3_mult_vec(_hkl2real, &pos);
	return pos;
}

AtomGroupPtr Crystal::getAtomsInBox(vec3 target, double tolx,
                                    double toly, double tolz)
{
	AtomGroupPtr atoms = AtomGroupPtr(new AtomGroup());

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
			
			atoms->addAtom(anAtom);
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

void Crystal::silentConcludeRefinement()
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();

	int num = _cycleNum;
	_silent = true;
	concludeRefinement(num, data);
	_silent = false;
}

void Crystal::wrapUpRefinement()
{
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();
	
	_lastRWork = _rWork;
	_cycleNum++;
	concludeRefinement(_cycleNum, data);
	saveState();

	if (_rWork <= _bestRWork)
	{
		_bestRWork = _rWork;
		_bestState = stateCount() - 1;
	}

	Options::getRuntimeOptions()->agreementSummary();
}

void Crystal::postRestoreState()
{
	Parser::postRestoreState();
	OptionsPtr options = Options::getRuntimeOptions();
	DiffractionPtr data = options->getActiveData();

	wrapUpRefinement();
}

void Crystal::openInCoot()
{
	std::string command = "coot " + _lastMtz + " " + _lastEnsemblePDB
	+ " " + _lastAveragePDB + "\n";
	
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
		Options::setNSamples(NULL, -1);
	}

	if (_sampleNum < 0) 
	{
		_sampleNum = 120;
	}

	double totalPoints = _sampleNum - 1;
	
	if (totalPoints < 0)
	{
		totalPoints = 0;
	}
	
	return totalPoints;
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
	
	for (int i = 0; i < _difft->nn(); i++)
	{
		vec3 frac = _difft->fracFromElement(i);
		
		if (frac.x > _spaceGroup->mapasu_zero[0] ||
		    frac.y > _spaceGroup->mapasu_zero[1] ||
		    frac.z > _spaceGroup->mapasu_zero[2]) 
		{
			continue;
		}
		
		Atom *atom = _bucket->nearbyAtom(i);
		double density = _difft->getReal(i);
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
	
	double sum_atoms = sum_all - sum_solvent;
	double nearAtom = sum_atoms / sum_all * 100;
	double inSolvent = sum_solvent / sum_all * 100;
	
	double cNearAtom = (nearAtom / num_atom);
	double cSolvent = (inSolvent / num_solvent);
	double cAll = cNearAtom + cSolvent;
	
	double backbone = sum_back / sum_atoms * 100;
	double sidechain = sum_side / sum_atoms * 100;
	double hetatm = sum_hetatm / sum_atoms * 100;

	/* Adjusted for volume of regions */
	double cBackbone = sum_back / num_back;
	double cSidechain = sum_side / num_side;
	double cHetatm = sum_hetatm / num_hetatm;
	
	if (cHetatm != cHetatm)
	{
		cHetatm = 0;
	}
	
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

double Crystal::getProbeRadius()
{
	if (Options::getProbeRadius() < 0 && _probeRadius <= 0)
	{
		_probeRadius = 0.4;
		Options::setProbeRadius(NULL, _probeRadius);
	}

	if (Options::getProbeRadius() >= 0)
	{
		_probeRadius = Options::getProbeRadius();
		Options::setProbeRadius(NULL, _probeRadius);
	}
	
	return _probeRadius;

}

double Crystal::getRealBFactor()
{
	if (Options::getGlobalBFactor() >= 0)
	{
		_realBFactor = Options::getGlobalBFactor();
		std::cout << "Picked up new B factor of " << _realBFactor <<
		"." <<  std::endl;
		Options::resetGlobalBFactor();
	}
	
	if (_realBFactor < 0)
	{
		_realBFactor = 0;
	}
	
	return _realBFactor;
}

void Crystal::updateLargestNum(AtomPtr atom)
{
	if (atom->getAtomNum() > _largestNum)
	{
		_largestNum = atom->getAtomNum();
	}
}

void Crystal::whack()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		MoleculePtr mol = molecule(i);
		
		if (!mol->isPolymer())
		{
			continue;
		}
		
		ToPolymerPtr(mol)->whack();
	}
}

void Crystal::chelate()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		MoleculePtr mole = molecule(i);

		if (mole->isPolymer())
		{
			continue;
		}
		
		mole->chelate("ZN", 15.);
		mole->chelate("CA", 15.);
		
		if (!mole->isWaterNetwork())
		{
			continue;
		}
	}
}

void Crystal::scaleAnchorBs(double ratio)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		MoleculePtr mole = molecule(i);

		if (!mole->isPolymer())
		{
			continue;
		}
		
		PolymerPtr pol = ToPolymerPtr(mole);
		AnchorPtr anch = pol->getAnchorModel();
		double b = anch->getBFactor();
		b *= ratio;
		anch->setBFactor(b);
		
		pol->propagateChange();
	}
}

double Crystal::averageBFactor()
{
	double ave = 0;
	double count = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		MoleculePtr mole = molecule(i);

		if (!mole->isPolymer())
		{
			continue;
		}
		
		ave += ToPolymerPtr(mole)->getAverageBFactor();
		count++;
	}
	
	return ave/count;
}

void Crystal::savePositions()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->saveAtomPositions();
	}
}

void Crystal::refitToSavedPositions()
{
	std::cout << "Refitting positions." << std::endl;
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->isPolymer())
		{
			PolymerPtr pol = ToPolymerPtr(molecule(i));
			for (int j = 0; j < 1; j++)
			{
				pol->getAnchorModel()->propagateChange(-1, true);
				pol->refineAnchorPosition(shared_from_this());
				pol->refine(shared_from_this(), RefinementSavedPos);
			}
		}
	}
}

void Crystal::removeAtom(AtomPtr atom)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->removeAtom(atom);
	}
}

void Crystal::refreshAnchors()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->isPolymer())
		{
			PolymerPtr pol = ToPolymerPtr(molecule(i));

			if (pol->getAnchorModel())
			{
				pol->getAnchorModel()->forceRefresh();
			}
		}
	}
}

void Crystal::addPDBContents(std::string pdbName)
{
	PDBReader pdb = PDBReader();
	pdb.setFilename(pdbName);
	CrystalPtr crystal = pdb.getCrystal();
	crystal->setAnchors();
	crystal->tieAtomsUp();
	crystal->hydrogenateContents();

	for (int i = 0; i < crystal->moleculeCount(); i++)
	{
		MoleculePtr mol = crystal->molecule(i);
		addMolecule(mol);
		vec3 centroid = mol->initialCentroid();
		Options::getRuntimeOptions()->focusOnPosition(centroid);
		
		if (!mol->isPolymer())
		{
			continue;
		}
		
		PolymerPtr pol = ToPolymerPtr(mol);

		MotionPtr overall = getOverallMotion();
		
		if (overall)
		{
			overall->addToPolymer(pol);
		}
		
		for (int j = 0; j < 5; j++)
		{
			pol->refine(shared_from_this(), RefinementModelPos);
		}

		pol->getAnchorModel()->atLeastOneMotion();
		_motions.back()->refine();
	}
}

void Crystal::resetMotions()
{
	for (int i = 0; i < _motions.size(); i++)
	{
		_motions[i]->reset();
	}
	
	refreshPositions();
}

void Crystal::addMotion(MotionPtr mot, PolymerPtr origPol)
{
	_motions.push_back(mot);
	
	if (!origPol)
	{
		return;
	}
	
	char first = origPol->getChainID()[0];

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer() || molecule(i) == origPol)
		{
			continue;
		}
		
		if (molecule(i)->getChainID()[0] == first)
		{
			mot->addToPolymer(ToPolymerPtr(molecule(i)));
		}
	}
}
