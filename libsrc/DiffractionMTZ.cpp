//
//  DiffractionMTZ.cpp
//  vagabond
//
//  Created by Helen Ginn on 22/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "DiffractionMTZ.h"
#include <iostream>
#include "Shouter.h"
#include "../libccp4/csymlib.h"
#include <vector>
#include "Options.h"
#include "FileReader.h"

void getCol(std::vector<std::string> names, CMtz::MTZ *mtz,
            CMtz::MTZCOL **column)
{

	for (size_t i = 0; i < names.size(); i++)
	{
		*column = MtzColLookup(mtz, names[i].c_str());

		if (*column)
		{
			std::cout << "Found column: " << names[i] << std::endl;
			break;
		}
	}
}

LabelChoice DiffractionMtz::prepareChoice(CMtz::MTZ *mtz)
{
	LabelChoice choice;

	int nxtals = MtzNxtal(mtz);

	for (size_t k = 0; k < nxtals; k++)
	{
		CMtz::MTZXTAL *xtal = MtzIxtal(mtz, k);
		int nsets = MtzNsetsInXtal(xtal);

		for (size_t j = 0; j < nsets; j++)
		{
			CMtz::MTZSET *set = MtzIsetInXtal(xtal, j);
			int ncol = MtzNcolsInSet(set);

			for (size_t i = 0; i < ncol; i++)
			{
				CMtz::MTZCOL *col = MtzIcolInSet(set, i);
				std::string label = col->label;
				std::string type = col->type;

				choice.availabels.push_back(label);
				choice.types.push_back(type);
			}
		}
	}
	
	return choice;
}

void DiffractionMtz::load()
{
	syminfoCheck();

	/* Assume that the filename is an MTZ file! */

	CMtz::MTZ *mtz = CMtz::MtzGet(_filename.c_str(), 0);
	//	int spgNum = MtzSpacegroupNumber(mtz);
	
	if (!file_exists(_filename))
	{
		Shouter *shout = new Shouter("MTZ file " + _filename 
		                             + " does not exist!");
		throw shout;
	}

	if (mtz == NULL)
	return;

	float *refldata;
	refldata= (float *) CCP4::ccp4_utils_malloc((mtz->ncol_read + 1)
	                                            * mtz->nref_filein
	* sizeof(float));

	memset(refldata, '\0', (mtz->ncol_read + 1) *
	       mtz->nref_filein * sizeof(float));

	int xtalCount = CMtz::MtzNxtal(mtz);
	std::cout << "Total crystals in " << _filename << ": " << 
	xtalCount << std::endl;
	
	CMtz::MTZXTAL **xtals = CMtz::MtzXtals(mtz);
	
	char *spgname = mtz->mtzsymm.spcgrpname;
	CSym::CCP4SPG *spg = CSym::ccp4spg_load_by_ccp4_spgname(spgname);

	if (xtalCount == 0)
	{
		std::cout << "Warning! No crystals in MTZ file?" << std::endl;
	}
	else
	{
		for (int i = 0; i < xtalCount; i++)
		{
			std::cout << "Crystal " << i + 1 << " unit cell dims: ";
			std::cout << "a = " << xtals[i]->cell[0] << " Å, ";
			std::cout << "b = " << xtals[i]->cell[1] << " Å, ";
			std::cout << "c = " << xtals[i]->cell[2] << " Å, ";
			std::cout << "alpha = " << xtals[i]->cell[3] << "°, ";
			std::cout << "beta = " << xtals[i]->cell[4] << "°, ";
			std::cout << "gamma = " << xtals[i]->cell[5] << "°." << std::endl;
		}
	}

	CMtz::MtzRrefl(mtz->filein, mtz->ncol_read * mtz->nref_filein, refldata);

	CMtz::MTZCOL *col_f = NULL;
	CMtz::MTZCOL *col_fwt = NULL;
	CMtz::MTZCOL *col_sigf = NULL;
	CMtz::MTZCOL *col_rfree = NULL;
	CMtz::MTZCOL *col_fpart = NULL;
	CMtz::MTZCOL *col_phase = NULL;
	CMtz::MTZCOL *col_phipart = NULL;

	std::vector<std::string> ampNames;
	std::string optAmp = Options::getLabF();
	if (optAmp.length())
	{
		ampNames.push_back(optAmp);
	}
	ampNames.push_back("F");
	ampNames.push_back("F-obs");
	ampNames.push_back("F-obs-filtered");
	ampNames.push_back("FP");

	std::vector<std::string> phNames, fwtNames;
	std::string optPh = Options::getLabPhase();
	std::string optFwt = Options::getLabFWT();
	if (optPh.length())
	{
		phNames.push_back(optPh);
	}
	phNames.push_back("PHWT");
	getCol(phNames, mtz, &col_phase);

	if (optFwt.length())
	{
		fwtNames.push_back(optFwt);
	}
	fwtNames.push_back("FWT");
	getCol(fwtNames, mtz, &col_fwt);
	
	std::vector<std::string> fParts, phiParts;
	fParts.push_back("Fpart");
	fParts.push_back("FPART");

	getCol(ampNames, mtz, &col_f);
	LabelChoice choice = prepareChoice(mtz);

	if (!col_f)
	{
		choice.original = "FP";
		choice.wanted = "observed amplitudes";
		Shouter *shout;
		shout = new Shouter("I could not find your amplitude column in\n"
		              + _filename + " - please label as F or FP.");
		shout->setChoice(choice);
		throw shout;
	}

	std::vector<std::string> errNames;
	
	std::string optSigF = Options::getLabSigF();
	if (optSigF.length())
	{
		errNames.push_back(optSigF);
	}
	
	errNames.push_back("SIGF");
	errNames.push_back("SIGFP");
	errNames.push_back("SIGF-obs");
	errNames.push_back("SIGF-obs-filtered");

	getCol(errNames, mtz, &col_sigf);

	if (!col_sigf)
	{
		choice.original = "SIGFP";
		choice.wanted = "sigma values on amplitudes";
		Shouter *shout;
		shout = new Shouter("I could not find your sigma/error column in\n"
		              + _filename + " - please label as SIGF or SIGFP.");
		shout->setChoice(choice);
		throw shout;
	}

	std::vector<std::string> rFreeNames;

	std::string optFree = Options::getLabFree();
	if (optFree.length())
	{
		rFreeNames.push_back(optFree);
	}

	rFreeNames.push_back("RFREE");
	rFreeNames.push_back("FREER");
	rFreeNames.push_back("FREE");
	rFreeNames.push_back("FreeR_flag");
	rFreeNames.push_back("R-free-flags");

	getCol(rFreeNames, mtz, &col_rfree);

	if (!col_rfree && _needsFree)
	{
		choice.original = "FREE";
		choice.wanted = "free set values";
		Shouter *shout;
		shout = new Shouter("I could not find your R-free-flag column in"
		              + _filename + " - please label as FREE or FreeRflag.");
		shout->setChoice(choice);
		throw shout;
	}

	CMtz::MTZCOL *col_h = MtzColLookup(mtz, "H");
	CMtz::MTZCOL *col_k = MtzColLookup(mtz, "K");
	CMtz::MTZCOL *col_l = MtzColLookup(mtz, "L");
	
	CMtz::MTZCOL *sorted[] = {col_h, col_k, col_l};
	
	CMtz::MtzSetSortOrder(mtz, sorted);

	if (!col_h || !col_k || !col_l)
	{
		Shouter *shout;
		shout = new Shouter("I could not find some or all of your\n"\
		                    "HKL indices columns " + _filename + 
		                    " - please\nlabel as H, K and L.");
		throw shout;
	}

	MtzResLimits(mtz, &_minRes, &_maxRes);
	_minRes = 1.0 / sqrt(_minRes);
	_maxRes = 1.0 / sqrt(_maxRes);
	
	if (_limit >= 0)
	{
		_maxRes = _limit;
	}
	
	bool uniform_cutoff = (_limit < 0);

	int indexLimitH = 0; int indexLimitK = 0; int indexLimitL = 0;

	for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
	{
		float *ptr = &refldata[i];
		int h = abs(ptr[col_h->source - 1]);
		int k = abs(ptr[col_k->source - 1]);
		int l = abs(ptr[col_l->source - 1]);
		
		indexLimitH = std::max(indexLimitH, h);
		indexLimitK = std::max(indexLimitK, k);
		indexLimitL = std::max(indexLimitL, l);
	}
		
		double hres = xtals[0]->cell[0] / indexLimitH;
		if (!uniform_cutoff && 1 / hres < _limit)
		{
			indexLimitH = xtals[0]->cell[0] / _limit + 1;
		}
		double kres = xtals[0]->cell[1] / indexLimitK;
		if (!uniform_cutoff && 1 / kres < _limit)
		{
			indexLimitK = xtals[0]->cell[1] / _limit + 1;
		}
		double lres = xtals[0]->cell[2] / indexLimitL;
		if (!uniform_cutoff && 1 / lres < _limit)
		{
			indexLimitL = xtals[0]->cell[2] / _limit + 1;
		}

		int largest = std::max(indexLimitH, indexLimitK);
		largest = std::max(largest, indexLimitL);

		largest = 2 * largest + 1;

	if (largest == 0)
	{
		Shouter *shout;
		shout = new Shouter("Problem determining resolution limits.\n"\
		              "Do you have a unit cell and some reflections?");
	}
	
	std::vector<double> unitCell;
	for (int i = 0; i < 6; i++)
	{
		unitCell.push_back(xtals[0]->cell[i]);
	}

	if (uniform_cutoff)
	{
		_fft = VagFFTPtr(new VagFFT(largest, largest, largest, 0, 1));
	}
	else
	{
		indexLimitH *= 2;
		indexLimitK *= 2;
		indexLimitL *= 2;

		_fft = VagFFTPtr(new VagFFT(indexLimitH, indexLimitK, 
		                            indexLimitL, 0, 1));
	}
	_fft->setUnitCell(unitCell);
	_fft->setSpaceGroup(spg);
	_fft->multiplyAll(nan(" "));
	_fft->copyToScratch(0);
	_fft->multiplyAll(nan(" "));

	int count = 0;
	int maskCount = 0;

	phiParts.push_back("PHIpart");
	phiParts.push_back("PHIPART");
	getCol(fParts, mtz, &col_fpart);
	getCol(phiParts, mtz, &col_phipart);
	
	bool addPartial = (col_fpart != NULL && col_phipart != NULL);
	
	if (Options::usePartial() && !addPartial)
	{
		std::string warn = "Asked to scale partial data set, but no\n"
		"partial data set found in FPART/PHIPART columns.";
		warn_user(warn);
	}
	
	if (!Options::usePartial())
	{
		addPartial = false;
	}
	
	if (addPartial)
	{
		std::cout << std::endl << "Found partial F/PHI columns; loading"
		" partial structure." << std::endl;
		_partial = VagFFTPtr(new VagFFT(largest, largest, largest));
		_partial->multiplyAll(nan(" "));
	}

	if (col_phase != NULL)
	{
		std::cout << std::endl << "Found original phases from MTZ." 
		<< std::endl;
		_original = VagFFTPtr(new VagFFT(largest, largest, largest));
		_original->multiplyAll(0);
		_original->setUnitCell(unitCell);
		_original->setSpaceGroup(spg);
		_original->setStatus(FFTReciprocalSpace);
	}

	for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
	{
		float *adata = &refldata[i];

		int _h = adata[col_h->source - 1];
		int _k = adata[col_k->source - 1];
		int _l = adata[col_l->source - 1];
		int h, k, l;
		
		CSym::ccp4spg_put_in_asu(spg, _h, _k, _l, &h, &k, &l);
		
		if (!_fft->withinBounds(h, k, l))
		{
			continue;
		}
		
		float amplitude = adata[col_f->source - 1];
		float phase = 0;
		float fwt = 0;

		
		if (col_phase != NULL)
		{
			phase = adata[col_phase->source - 1];
		}
		
		if (col_fwt != NULL)
		{
			fwt = adata[col_fwt->source - 1];
		}

		float sigma = amplitude / 10;
		
		if (col_sigf != NULL)
		{
			sigma = adata[col_sigf->source - 1];
		}

		float flag = 1;

		if (col_rfree)
		{
			flag = adata[col_rfree->source - 1];
		}

		if (flag != flag) flag = 1;

		if (flag == 0 && amplitude != amplitude) flag = 1;

		float mask = flag;

		long element = _fft->element(h, k, l);

		_fft->setComponent(element, 0, amplitude);
		_fft->setComponent(element, 1, sigma);
		_fft->setScratchComponent(element, 0, 0, mask + 0.1);

		count++;

		if (mask < 0.5)
		{
			maskCount++;
		}

		if (addPartial)
		{
			float fpart = adata[col_fpart->source - 1];
			float phipart = adata[col_phipart->source - 1];
			phipart = deg2rad(phipart);

			float x = fpart * cos(phipart);
			float y = fpart * sin(phipart);

			_partial->setComponent(element, 0, x);
			_partial->setComponent(element, 1, y);
		}

		if (col_phase != NULL)
		{
			double ph = deg2rad(phase);
			float x = fwt * cos(ph);
			float y = fwt * sin(ph);
			_original->setComponent(element, 0, x);
			_original->setComponent(element, 1, y);
		}
	}

	std::cout << "Loaded " << count << " reflections into"\
	" memory from " << _filename << "." << std::endl;
	std::cout << "Counted " << maskCount << " free reflections." << std::endl << std::endl;

	if (_fft)
	{
		_fft->applySymmetry(true, -1, false);
	}

	free(refldata);
	MtzFree(mtz);
	
	return;

	/* Apply all symmetry relations */ 
	vec3 nLimits = getNLimits(_original, _original);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				int _h, _k, _l;
				int sym = CSym::ccp4spg_put_in_asu(spg, i, j, k, 
				                                   &_h, &_k, &_l);

				long index = _original->element(i, j, k);
				long asuIndex = _original->element(_h, _k, _l);

				double x = _original->getReal(asuIndex);
				double y = _original->getImag(asuIndex);

				/* establish amp & phase */
				double amp = sqrt(x * x + y * y);
				double myPhase = atan2(y, x);

				if (sym % 2 == 0)
				{
					myPhase *= -1;
				}

				int symop = (sym - 1) / 2;

				/* calculate phase shift for symmetry operator */
				float *trn = spg->symop[symop].trn;

				/* rotation */
				double shift = (float)i * trn[0];
				shift += (float)j * trn[1];
				shift += (float)k * trn[2];

				shift = fmod(shift, 1.);

				/*  apply shift when filling in other sym units */
				double newPhase = myPhase + shift * 2 * M_PI;

				x = amp * cos(newPhase);
				y = amp * sin(newPhase);

				if (_original)
				{
					_original->setComponent(index, 0, x);
					_original->setComponent(index, 1, y);
				}
			}
		}
	}
}

void DiffractionMtz::syminfoCheck()
{
	CSym::CCP4SPG *spg = CSym::ccp4spg_load_by_ccp4_num(1);
	if (!spg)
	{
		shout_at_user("I can't seem to load space group P1!\n" \
		              "If CCP4 is installed, it may need sourcing.\n"\
		              "Please make sure the value of the\n"\
		              "environment variable $SYMINFO is set correctly.\n"\
		              "Something like (bash):\n\texport SYMINFO=/path/to/ccp4-vX.X.X/lib/data/syminfo.lib\n"\
		              "Something like (csh):\n\tsetenv SYMINFO /path/to/ccp4-vX.X.X/lib/data/syminfo.lib\n"\
		              "\t(... make sure you correct these paths!)");

	}
}
