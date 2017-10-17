//
//  Crystal.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Crystal.h"
#include "fftw3d.h"
#include "vec3.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <time.h>
#include "BucketUniform.h"
#include "Shouter.h"
#include "Diffraction.h"
#include "Polymer.h"
#include "CSV.h"
#include "FileReader.h"

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
	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tieAtomsUp();
	}
}

void Crystal::addMolecule(MoleculePtr molecule)
{
	if (molecule->getChainID().length() <= 0)
	{
		shout_at_helen("Monomer chain ID is missing while trying\n"\
					   "to interpret PDB file.");
	}
	
	_molecules[molecule->getChainID()] = molecule;
}

void Crystal::setReal2HKL(mat3x3 mat)
{
	_real2frac = mat;
}

void Crystal::setHKL2Real(mat3x3 mat)
{
	_hkl2real = mat;
}

void Crystal::realSpaceClutter()
{
	if (!_fft)
	{
		_fft = FFTPtr(new FFT());
		_difft = FFTPtr(new FFT());

		vec3 uc_dims = empty_vec3();
		vec3 fft_dims = empty_vec3();
		uc_dims.x = mat3x3_length(_hkl2real, 0) / PROTEIN_SAMPLING;
		uc_dims.y = mat3x3_length(_hkl2real, 1) / PROTEIN_SAMPLING;
		uc_dims.z = mat3x3_length(_hkl2real, 2) / PROTEIN_SAMPLING;

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


	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->addToMap(_fft, _real2frac);
	}

//	BucketPtr bucket = BucketPtr(new BucketUniform());
//	bucket->addSolvent(fft);

	_fft->createFFTWplan(8);
	_difft->createFFTWplan(8);
}

void Crystal::writeCalcMillersToFile(DiffractionPtr data, std::string prefix)
{
	if (!_fft)
	{
		shout_at_user("There is likely a bug. Cannot write\n"\
					  "calculated Miller list until it has\n"\
					  "been generated!");
	}

	_fft->setAll(0);
	realSpaceClutter();
	fourierTransform(1);
	scaleToDiffraction(data);
	double nLimit = _fft->nx;
	nLimit = nLimit - ((int)nLimit % 2);
	nLimit /= 2;

	double dStar = 1 / _maxResolution;

	mat3x3 transpose = mat3x3_transpose(_hkl2real);

	double aLength = mat3x3_length(transpose, 0);
	double bLength = mat3x3_length(transpose, 1);
	double cLength = mat3x3_length(transpose, 2);

	double aLimit = aLength * dStar;
	double bLimit = bLength * dStar;
	double cLimit = cLength * dStar;

	/* For writing MTZ files */

	int columns = 10;

	float cell[6], wavelength;
	float *fdata = new float[columns];

	/* variables for symmetry */
	CSym::CCP4SPG *mtzspg = _spaceGroup;
	float rsm[192][4][4];
	char ltypex[2];

	/* variables for MTZ data structure */
	CMtz::MTZ *mtzout;
	CMtz::MTZXTAL *xtal;
	CMtz::MTZSET *set;
	CMtz::MTZCOL *colout[11];

	cell[0] = _unitCell[0];
	cell[1] = _unitCell[1];
	cell[2] = _unitCell[2];
	cell[3] = _unitCell[3];
	cell[4] = _unitCell[4];
	cell[5] = _unitCell[5];
	wavelength = 1.00; // fixme

	std::string outputFile = prefix + "_" + _filename + "_vbond.mtz";

	mtzout = CMtz::MtzMalloc(0, 0);
	ccp4_lwtitl(mtzout, "Written from Helen's XFEL tasks ", 0);
	mtzout->refs_in_memory = 0;
	mtzout->fileout = CMtz::MtzOpenForWrite(outputFile.c_str());

	// then add symm headers...
	for (int i = 0; i < mtzspg->nsymop; ++i)
		CCP4::rotandtrn_to_mat4(rsm[i], mtzspg->symop[i]);
	strncpy(ltypex, mtzspg->symbol_old, 1);
	ccp4_lwsymm(mtzout, mtzspg->nsymop, mtzspg->nsymop_prim, rsm, ltypex,
				mtzspg->spg_ccp4_num, mtzspg->symbol_old, mtzspg->point_group);

	// then add xtals, datasets, cols
	xtal = MtzAddXtal(mtzout, "vagabond_crystal", "vagabond_project", cell);
	set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
	colout[0] = MtzAddColumn(mtzout, set, "H", "H");
	colout[1] = MtzAddColumn(mtzout, set, "K", "H");
	colout[2] = MtzAddColumn(mtzout, set, "L", "H");
	colout[3] = MtzAddColumn(mtzout, set, "FREE", "R");
	colout[4] = MtzAddColumn(mtzout, set, "FP", "F");
	colout[5] = MtzAddColumn(mtzout, set, "FC", "F");
	colout[6] = MtzAddColumn(mtzout, set, "FWT", "F");
	colout[7] = MtzAddColumn(mtzout, set, "PHWT", "P");
	colout[8] = MtzAddColumn(mtzout, set, "DELFWT", "F");
	colout[9] = MtzAddColumn(mtzout, set, "PHDELWT", "P");

	int num = 0;
	FFTPtr fftData = data->getFFT();

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
			{
				bool asu = CSym::ccp4spg_is_in_asu(_spaceGroup, i, j, k);

				if (!asu)
				{
					continue;
				}

				vec3 pos = make_vec3(i, j, k);
				mat3x3_mult_vec(_real2frac, &pos);

				double phase = _fft->getPhase(i, j, k);

				double intensity = _fft->getIntensity(i, j, k);
				double calcAmp = sqrt(intensity);

				double foInt = fftData->getIntensity(i, j, k);
				double foAmp = sqrt(foInt);

				int free = fftData->getMask(i, j, k);

				if (vec3_length(pos) > dStar)
				{
					continue;
				}

				/* MTZ file stuff */

				fdata[0] = i;
				fdata[1] = j;
				fdata[2] = k;
				fdata[3] = free;
				fdata[4] = foAmp;
				fdata[5] = calcAmp;
				fdata[6] = 2 * foAmp - calcAmp;
				fdata[7] = phase;
				fdata[8] = foAmp - calcAmp;
				fdata[9] = phase;

				num++;
				ccp4_lwrefl(mtzout, fdata, colout, columns, num);
			}
		}
	}

	MtzPut(mtzout, " ");
	MtzFree(mtzout);

	std::cout << "Written pha files from crystal." << std::endl;
}

double Crystal::valueWithDiffraction(DiffractionPtr data, two_dataset_op op,
									 bool verbose, double lowRes, double highRes)
{
	if (!_fft || !_fft->nn)
	{
		realSpaceClutter();
		scaleToDiffraction(data);
	}

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
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(_real2frac, &ijk);
				double length = vec3_length(ijk);

				if (length < minRes || length > maxRes)
				{
					continue;
				}

				double amp1 = sqrt(fftData->getIntensity(i, j, k));
				double amp2 = sqrt(_fft->getIntensity(i, j, k));

				int isFree = fftData->getMask(i, j, k) == MaskFree;

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
		csv->writeToFile("correlplot.csv");

		std::map<std::string, std::string> plotMap;
		plotMap["filename"] = "correlplot";
		plotMap["xHeader0"] = "fo";
		plotMap["yHeader0"] = "fc";
		plotMap["colour0"] = "black";

		plotMap["xTitle0"] = "Fo amplitude";
		plotMap["yTitle0"] = "Fc amplitude";
		plotMap["style0"] = "scatter";
		csv->plotPNG(plotMap);
	}

	double working = (*op)(set1, set2);

	if (verbose)
	{
		double ccWork = correlation(set1, set2);
		double ccFree = correlation(free1, free2);
		double free = (*op)(free1, free2);
		double diff = free - working;

		std::cout << "CCwork/CCfree: " << ccWork * 100 << ", " << ccFree * 100
		<< " %." << std::endl;

		std::cout << "R values: " << std::setprecision(4)
		<< working * 100;
		std::cout << ", " << free * 100 << " % (diff: " << diff * 100 << " %)"<<  std::endl;
	}

	return working;
}

void Crystal::applyScaleFactor(double scale, double lowRes, double highRes)
{
	double nLimit = _fft->nx;
	nLimit /= 2;
	std::vector<double> set1, set2, free1, free2;

	double minRes = (lowRes == 0 ? 0 : 1 / lowRes);
	double maxRes = (highRes == 0 ? FLT_MAX : 1 / highRes);

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
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

void Crystal::scaleToDiffraction(DiffractionPtr data)
{
	std::vector<double> bins;
	generateResolutionBins(0, _maxResolution, 20, &bins);

	for (int i = 0; i < bins.size() - 1; i++)
	{
		double ratio = valueWithDiffraction(data, &scale_factor_by_sum, false,
											bins[i], bins[i + 1]);
		double scale = 1 / ratio;
		applyScaleFactor(scale, bins[i], bins[i + 1]);
	}

}

double Crystal::rFactorWithDiffraction(DiffractionPtr data, bool verbose)
{
	std::cout << "*******************************" << std::endl;

	double rFactor = valueWithDiffraction(data, &r_factor, verbose);

	std::cout << "*******************************" << std::endl;

	return rFactor;
}

void Crystal::getDataInformation(DiffractionPtr data, double partsFo,
								 double partsFc)
{
	if (!_fft || !_fft->nn)
	{
		realSpaceClutter();
	}

	fourierTransform(1);

	scaleToDiffraction(data);
	rFactorWithDiffraction(data, true);

	std::cout << "Combining model with data information." << std::endl;

	FFTPtr fftData = data->getFFT();
	double nLimit = std::min(fftData->nx, _fft->nx);
	nLimit /= 2;
//	nLimit += 1;
	std::vector<double> set1, set2;

	/* symmetry issues */
	for (int i = -nLimit; i < nLimit; i++)
	{
		for (int j = -nLimit; j < nLimit; j++)
		{
			for (int k = 0; k < nLimit; k++)
			{
				double amp = sqrt(fftData->getIntensity(i, j, k));
				bool isRfree = fftData->getMask(i, j, k) == MaskFree;

				if (amp != amp || isRfree)
				{
					continue;
				}

				vec2 complex;
				long index = _fft->element(i, j, k);
				complex.x = _fft->getReal(index);
				complex.y = _fft->getImaginary(index);
				double old_amp = sqrt(complex.x * complex.x +
									  complex.y * complex.y);

				double new_amp = partsFo * amp - partsFc * old_amp;
				new_amp /= old_amp;

				double diff_scale = amp - old_amp;
				diff_scale /= old_amp;

				vec2 diff_complex = complex;

				complex.x *= new_amp;
				complex.y *= new_amp;

				diff_complex.x *= diff_scale;
				diff_complex.y *= diff_scale;

				_fft->setElement(index, complex.x, complex.y);
				_difft->setElement(index, diff_complex.x, diff_complex.y);
			}
		}
	}

	fourierTransform(-1);
	_difft->fft(-1);
}

void Crystal::tiedUpScattering()
{
	double tied = 0;
	double total = 0;

	for (int i = 0; i < moleculeCount(); i++)
	{
		molecule(i)->tiedUpScattering(&tied, &total);
	}

	std::cout << "Tied up " << tied << " electrons out of " << total << " (";
	std::cout << 100. * sqrt(tied / total) << "%)." << std::endl;
}

void Crystal::setAnchors()
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() == "Polymer")
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));
			polymer->setAnchor(_anchorResidue);
		}
	}
}


void Crystal::changeAnchors(int newAnchor)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() == "Polymer")
		{
			PolymerPtr polymer = ToPolymerPtr(molecule(i));

			polymer->changeAnchor(newAnchor);
		}
	}
}

Crystal::Crystal()
{
	_firstScale = -1;
	_maxResolution = HARD_CODED_RESOLUTION;
	_anchorResidue = 34;
	_overallB = 16;
}

void Crystal::applySymOps()
{
	if (_spaceGroup->spg_num == 1)
	{
		return;
	}

	std::cout << "Applying space group " << _spaceGroup->symbol_xHM << std::endl;

//	return;
	_fft->applySymmetry(_spaceGroup, false);
	_fft->applySymmetry(_spaceGroup, true);
}

void Crystal::fourierTransform(int dir)
{
	_fft->fft(dir);

	if (dir == 1)
	{
		applySymOps();
	}
}

void Crystal::concludeRefinement(int cycleNum, DiffractionPtr data)
{
	for (int i = 0; i < moleculeCount(); i++)
	{
		if (molecule(i)->getClassName() != "Polymer")
		{
			continue;
		}

		realSpaceClutter();
		getDataInformation(data, 3, 2);
		std::string refineCount = "refine_" + i_to_str(cycleNum);
		PolymerPtr polymer = ToPolymerPtr(molecule(i));
		polymer->makePDB(refineCount + ".pdb");
		polymer->graph("graph_" + i_to_str(cycleNum));
		writeCalcMillersToFile(data, refineCount);
		polymer->differenceGraphs("diffgraph_" + i_to_str(cycleNum), shared_from_this());
	}
}
