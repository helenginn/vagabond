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
#include "../libccp4/cmtzlib.h"
#include <vector>
#include "fftw3d.h"
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
			break;
		}
	}
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

	float *adata = (float *) CCP4::ccp4_utils_malloc((mtz->ncol_read)
	                                                 * sizeof(float));
	
	int xtalCount = CMtz::MtzNxtal(mtz);
	std::cout << "Total crystals in " << _filename << ": " << 
	xtalCount << std::endl;
	
	CMtz::MTZXTAL **xtals = CMtz::MtzXtals(mtz);

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
	CMtz::MTZCOL *col_sigf = NULL;
	CMtz::MTZCOL *col_rfree = NULL;
	CMtz::MTZCOL *col_fpart = NULL;
	CMtz::MTZCOL *col_phipart = NULL;

	std::vector<std::string> ampNames;
	ampNames.push_back("F");
	ampNames.push_back("FP");
	
	std::vector<std::string> fParts, phiParts;
	fParts.push_back("Fpart");
	fParts.push_back("FPART");

	getCol(ampNames, mtz, &col_f);
	if (!col_f)
	{
		Shouter *shout;
		shout = new Shouter("I could not find your amplitude column in\n"
		              + _filename + " - please label as F or FP.");
		throw shout;
	}

	std::vector<std::string> errNames;
	errNames.push_back("SIGF");
	errNames.push_back("SIGFP");

	getCol(errNames, mtz, &col_sigf);
	if (!col_sigf && false)
	{
		warn_user("I could not find your sigma/error column in\n"
		          + _filename + " - please label as SIGF or SIGFP.\n"
		"I can do without and will keep going.");
	}


	std::vector<std::string> rFreeNames;
	rFreeNames.push_back("RFREE");
	rFreeNames.push_back("FREER");
	rFreeNames.push_back("FREE");
	rFreeNames.push_back("FreeR_flag");

	getCol(rFreeNames, mtz, &col_rfree);

	if (!col_rfree)
	{
		warn_user("I could not find your R-free-flag column in\n"
		          + _filename + " - please label as FREE or RFREE.\n"
		"I can do without and will keep going.");
	}

	CMtz::MTZCOL *col_h = MtzColLookup(mtz, "H");
	CMtz::MTZCOL *col_k = MtzColLookup(mtz, "K");
	CMtz::MTZCOL *col_l = MtzColLookup(mtz, "L");

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

	int largest = 0;

	int indexLimitH = std::max(fabs(col_h->min), fabs(col_h->max));
	int indexLimitK = std::max(fabs(col_k->min), fabs(col_k->max));
	int indexLimitL = std::max(fabs(col_l->min), fabs(col_l->max));

	largest = std::max(indexLimitH, indexLimitK);
	largest = std::max(largest, indexLimitL);

	largest *= 2;
	largest += 1;

	if (largest == 0)
	{
		Shouter *shout;
		shout = new Shouter("Problem determining resolution limits.\n"\
		              "Do you have a unit cell and some reflections?");
	}

	fft = FFTPtr(new FFT());
	fft->create(largest);
	fft->multiplyAll(nan(" "));
	fft->setupMask();

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
		_partial = FFTPtr(new FFT());
		_partial->create(largest);
		_partial->multiplyAll(nan(" "));
	}

	for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
	{
		memcpy(adata, &refldata[i], mtz->ncol_read * sizeof(float));

		int h = adata[col_h->source - 1];
		int k = adata[col_k->source - 1];
		int l = adata[col_l->source - 1];
		float amplitude = adata[col_f->source - 1];
		float sigma = adata[col_sigf->source - 1];
		float flag = 1;

		if (col_rfree)
		{
			flag = adata[col_rfree->source - 1];
		}

		if (flag != flag) flag = 1;
		
		if (flag == 0 && amplitude != amplitude) flag = 1;

		int mask = flag;

		long element = fft->element(h, k, l);

		fft->data[element][0] = amplitude;
		fft->data[element][1] = sigma;
		fft->mask[element] = mask;

		count++;

		if (mask == 0)
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
			
			_partial->data[element][0] = x;
			_partial->data[element][1] = y;
		}
	}
	
	free(adata);

	std::cout << "Loaded " << count << " reflections into"\
	" memory from " << _filename << "." << std::endl;
	std::cout << "Counted " << maskCount << " free reflections." << std::endl << std::endl;

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
