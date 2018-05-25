//
//  Crystal.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017-8 Helen Ginn. All rights reserved.
//

#include "Crystal.h"
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
#include "LocalCC.h"
#include "Atom.h"
#include "Kabsch.h"
#include "RefinementGridSearch.h"
#include "BoneDensity.h"
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
			if (maxRes >= 2.5)
			{
				sampling = maxRes / 5.;
			}
			if (maxRes <= 1.2)
			{
				sampling = maxRes / 3.;
			}
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

		fft_dims.x = largest; fft_dims.y = largest; fft_dims.z = largest;

		_fft->create(fft_dims.x, fft_dims.y, fft_dims.z);
		_fft->setupMask();

		_difft->create(fft_dims.x, fft_dims.y, fft_dims.z);

		double scaling = 1 / largest;

		_fft->setBasis(_hkl2real, scaling);
		_difft->setBasis(_hkl2real, scaling);
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

void Crystal::writeMillersToFile(DiffractionPtr data, std::string prefix)
{
	std::string outputFileOnly = prefix + "_" + _filename + "_vbond.mtz";
	getFFT()->writeReciprocalToFile(outputFileOnly, _maxResolution, _spaceGroup,
	                                _unitCell, _real2frac, data->getFFT());
	std::string outputFile = FileReader::addOutputDirectory(outputFileOnly);
	
	_lastMtz = outputFile;
	
	if (_bucket)
	{
		_bucket->writeMillersToFile(prefix, _maxResolution);	
	}
}

double Crystal::valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
                                     bool verbose, double lowRes, double highRes)
{
	FFTPtr fftData = data->getFFT();
	double nLimit = std::min(fftData->nx, _fft->nx);
	nLimit = nLimit - ((int)nLimit % 2);
	nLimit /= 2;

	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes == 0 ? 1 / _maxResolution : 1 / highRes);

	CSVPtr csv = CSVPtr(new CSV(2, "fo" , "fc"));

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
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

void Crystal::applyScaleFactor(double scale, double lowRes, double highRes)
{
	double nLimit = _fft->nx;
	nLimit /= 2;
	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes <= 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes <= 0 ? FLT_MAX : 1 / highRes);

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = -nLimit; k < nLimit; k++)
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

				if (real != real || imag != imag)
				{
					continue;
				}

				real *= scale;
				imag *= scale;

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
		_maxResolution = data->getMaxResolution();
		std::cout << std::setprecision(2);
		std::cout << "Using the resolution from " << data->getFilename()
		<< " of " << _maxResolution << " Å." << std::endl;
	}
	
	return _maxResolution;
}

/* bigger number is more detail */
double Crystal::getMaximumDStar(DiffractionPtr data)
{
	double maxRes = getMaxResolution(data);
	maxRes = 1 / maxRes;
	maxRes *= 1.4;

	return maxRes;
}

void Crystal::scaleToDiffraction(DiffractionPtr data)
{
	getMaxResolution(data);
	
	/* First, apply a scale factor to the entire range */
	double totalFc = totalToScale();
	double ratio = valueWithDiffraction(data, &scale_factor_by_sum, false,
	                                    0, _maxResolution);
	applyScaleFactor(totalFc / ratio, 0, 0);

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

void Crystal::scaleComponents(DiffractionPtr data)
{
	scaleToDiffraction(data);
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
	
	fourierTransform(1, data->getMaxResolution());
	scaleComponents(data);
	
	writeMillersToFile(data, prefix);

	double rFac = rFactorWithDiffraction(data, true);

	FFTPtr fftData = data->getFFT();
	double nLimit = std::min(fftData->nx, _fft->nx);
	nLimit /= 2;
	std::vector<double> set1, set2;

	double lowRes = Options::minRes();
	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = (1 / _maxResolution);
	
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
		std::cout << "There were " << bad << " bad reflections";
		std::cout << " out of " << _fft->nn << "!" << std::endl;
	}

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = -nLimit; k < nLimit; k++)
			{
				int _h, _k, _l;
				CSym::ccp4spg_put_in_asu(_spaceGroup, i, j, k, &_h, &_k, &_l);
				int isAbs = CSym::ccp4spg_is_sysabs(_spaceGroup, i, j, k);

				double amp = sqrt(fftData->getIntensity(_h, _k, _l));
				bool isRfree = (fftData->getMask(_h, _k, _l) == 0);
				long index = _fft->element(i, j, k);

				vec3 ijk = make_vec3(i, j, k);    
				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				if (length < minRes || length > maxRes
				    || (isRfree && amp == amp) || isAbs)    
				{	
					_fft->setElement(index, 0, 0);
					_difft->setElement(index, 0, 0);

					continue;
				}

				if (amp != amp)
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

	/* Back to real space */
	fourierTransform(-1);
	_difft->fft(-1);

	CSVPtr csv = CSVPtr(new CSV(3, "real_obs", "real_calc", "solvent"));
	std::vector<double> real_mixed, chosen_calc;
	int count = 0;
	
	for (int i = 0; i < _fft->nn; i += skip)
	{
		double obs = _fft->data[i][0];
		double calc = real_calcs[count];
		count++;
		if (calc <= 0) continue;

		double solvent = _bucket->isSolvent(i);

		real_mixed.push_back(obs);
		chosen_calc.push_back(calc);
		csv->addEntry(3, obs, calc, solvent);
	}
	
	double correl = correlation(chosen_calc, real_mixed);
	csv->writeToFile("real_space_cc.csv");
	printf("Real space correlation coefficient: %.3f\n", correl);
	
	_bucket->abandonCalculations();
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

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() == "Polymer")
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			if (!_anchorResidues.size())
			{
				polymer->findAnchorNearestCentroid();
			}
			else
			{
				polymer->setAnchor(_anchorResidues[0]);
			}
		}
	}
}

Crystal::Crystal()
{
	_cycleNum = 0;
	_lastRWork = FLT_MAX;
	_bestRWork = FLT_MAX;
	_sinceBestNum = 0;
	_correlPlotNum = 0;
	_tied = false;
	_spaceGroup = NULL;
	_spgNum = 0;
	_maxResolution = 0;
	_solvScale = 0.5;
	_solvBFac = 10;
	_unitCell.resize(6);
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

void Crystal::makePDBs(std::string suffix)
{
	std::vector<std::string> prefices; std::vector<PDBType> pdbTypes;
	prefices.push_back("ensemble_"); pdbTypes.push_back(PDBTypeEnsemble);
	prefices.push_back("average_"); pdbTypes.push_back(PDBTypeAverage);

	for (int i = 0; i < prefices.size(); i++)
	{
		std::string path;
		path = FileReader::addOutputDirectory(prefices[i] + suffix + ".pdb");
		
		if (i == 0)
		{
			_lastEnsemblePDB = path;
		}
		else
		{
			_lastAveragePDB = path;
		}
		
		std::ofstream file;
		file.open(path);

		for (int j = 0; j < moleculeCount(); j++)
		{
			CrystalPtr crystal = shared_from_this();
			file << molecule(j)->makePDB(pdbTypes[i], crystal); 
		}

		file.close();
	}
};

void Crystal::writeVagabondFile(int cycleNum)
{
	std::ofstream file;
	std::string filename = "refine_" + i_to_str(cycleNum) + ".vbond";
	std::string vbondFile = FileReader::addOutputDirectory(filename);
	file.open(vbondFile);
	writeToFile(file, 0);
	file.close();
}

double Crystal::concludeRefinement(int cycleNum, DiffractionPtr data)
{
	_cycleNum = cycleNum;
	std::cout << "*******************************" << std::endl;
	std::cout << "\tCycle " << cycleNum << std::endl;

	std::string refineCount = "refine_" + i_to_str(cycleNum);

	double rFac = getDataInformation(data, 2, 1, refineCount);
	makePDBs(refineCount);

	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() == "Polymer")
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			if (cycleNum > 0)
			{
				polymer->differenceGraphs("density_" + polymer->getChainID() +
				                          "_" + i_to_str(cycleNum), shared_from_this());
			}
			polymer->graph("chain_" + polymer->getChainID() +
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

	return rFac;
}

void Crystal::reconfigureUnitCell()
{
	PolymerPtr polymer = ToPolymerPtr(molecule(0));

	Kabsch kabsch;

	std::vector<vec3> xs, ys;

	for (int i = 0; i < polymer->atomCount(); i++)
	{
		if (!polymer->atom(i)->isBackbone())
		{
			continue;
		}

		vec3 initPos = polymer->atom(i)->getPDBPosition();
		vec3 nowPos = polymer->atom(i)->getAbsolutePosition();

		xs.push_back(nowPos);
		ys.push_back(initPos);
	}

	kabsch.setAtoms(xs, ys);
	kabsch.fixCentroids();
	kabsch.run();

	mat3x3 transform = kabsch.findFinalTransform();
	//mat3x3 invTrans = mat3x3_inverse(transform);

	std::cout << mat3x3_desc(transform) << std::endl;

	mat3x3 potential = mat3x3_mult_mat3x3(transform, _hkl2real);

	double vals[6];
	unit_cell_from_mat3x3(potential, vals);

	std::cout << "Potential new unit cell: " << std::endl;
	std::cout << "\t" << std::endl;

	for (int i = 0; i < 6; i++)
	{
		std::cout << vals[i] << " ";
	}

	std::cout << std::endl;
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
	ss << "CCwork/free: " << _ccWork << ", " << _ccFree << std::endl;
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

double Crystal::vsRefineBackboneToDensity(void *object)
{
	Parser *parser = static_cast<Parser *>(object);
	Crystal *crystal = dynamic_cast<Crystal *>(parser);

	crystal->backboneDensityAnalysis();
}

void Crystal::backboneDensityAnalysis()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (!molecule(i)->isPolymer())
		{
			continue;
		}
		
		ToPolymerPtr(molecule(i))->refineBackbone();
		
		continue;

		BoneDensity density;
		density.setCrystal(shared_from_this());
		PolymerPtr polymer = ToPolymerPtr(molecule(i));
		density.setPolymer(polymer);
		density.analyse();
		
		for (int j = 0; j < density.instructionCount(); j++)
		{
			BackboneInstruction inst = density.instruction(j);

			std::cout << "Refining from " << inst.startRes << " to " <<
			inst.endRes << " by ";
			std::cout << (inst.rType == RefinementFine ? "correlation."
			: "squeezing.") << std::endl;
			polymer->clearParams();

			switch (inst.rType)
			{
				case RefinementFine:
				polymer->addParamType(ParamOptionTorsion, 0.02);
				polymer->addParamType(ParamOptionKick, 0.010);
				polymer->addParamType(ParamOptionDampen, 0.005);
				polymer->addParamType(ParamOptionMagicAngles, 3);
				polymer->addParamType(ParamOptionNumBonds, 8);
				break;
				
				case RefinementModelRMSDZero:
				polymer->addParamType(ParamOptionTorsion, 0.04);
				polymer->addParamType(ParamOptionDampen, 0.005);
				polymer->addParamType(ParamOptionMagicAngles, 3);
				polymer->addParamType(ParamOptionNumBonds, 12);
				break;
				
				default:
				break;
			}
			
			polymer->refineRange(inst.startRes, inst.endRes,
			                     shared_from_this(), inst.rType);
		}
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

		ToPolymerPtr(molecule(i))->optimiseWholeMolecule(translation, rotation);
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
	addIntProperty("cycles_since_best", &_sinceBestNum);

	_spgNum = 0;
	if (_spaceGroup)
	{
		_spgNum = _spaceGroup->spg_num;
	}

	addIntProperty("spacegroup", &(_spgNum));

	for (int i = 0; i < moleculeCount(); i++)
	{
		addChild("molecule", molecule(i));
	}
	
	exposeFunction("recalculate_map", Crystal::vsConcludeRefinement);
	exposeFunction("restore_state", Crystal::vsRestoreState);
	exposeFunction("refine_backbone_to_density",
	               Crystal::vsRefineBackboneToDensity);
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
	_spaceGroup = CSym::ccp4spg_load_by_ccp4_num(_spgNum);
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

