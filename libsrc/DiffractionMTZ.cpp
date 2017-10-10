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

void getCol(std::vector<std::string> names, CMtz::MTZ *mtz,
			CMtz::MTZCOL **column)
{

	for (int i = 0; i < names.size(); i++)
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

	CMtz::MtzRrefl(mtz->filein, mtz->ncol_read * mtz->nref_filein, refldata);

	CMtz::MTZCOL *col_f = NULL;
	CMtz::MTZCOL *col_sigf = NULL;
	CMtz::MTZCOL *col_rfree = NULL;

	std::vector<std::string> ampNames;
	ampNames.push_back("F");
	ampNames.push_back("FP");

	getCol(ampNames, mtz, &col_f);
	if (!col_f)
	{
		shout_at_user("I could not find your amplitude column in\n"
					  + _filename + " - please label as F or FP.");
	}

	std::vector<std::string> errNames;
	errNames.push_back("SIGF");
	errNames.push_back("SIGFP");

	getCol(errNames, mtz, &col_sigf);
	if (!col_sigf)
	{
		warn_user("I could not find your sigma/error column in\n"
				  + _filename + " - please label as SIGF or SIGFP.\n"
				  "I can do without and will keep going.");
	}


	std::vector<std::string> rFreeNames;
	rFreeNames.push_back("RFREE");
	rFreeNames.push_back("FREE");

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
		shout_at_user("I could not find some or all of your\n"\
					  "HKL indices columns " + _filename + " - please\n"\
					  "label as H, K and L.");
	}

	float minRes, maxRes;
	MtzResLimits(mtz, &minRes, &maxRes);

	CMtz::MTZXTAL **xtals = MtzXtals(mtz);
	float *cell = (float *)malloc(sizeof(float) * 6);
	ccp4_lrcell(xtals[0], cell);

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
		shout_at_user("Problem determining resolution limits.\n"\
					  "Do you have a unit cell and some reflections?");
	}

	fft = FFTPtr(new FFT());
	fft->create(largest);
	fft->multiplyAll(nan(" "));
	fft->setupMask();

	int count = 0;
	int maskCount = 0;

	for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
	{
		memcpy(adata, &refldata[i], mtz->ncol_read * sizeof(float));

		int h = adata[col_h->source - 1];
		int k = adata[col_k->source - 1];
		int l = adata[col_l->source - 1];
		float amplitude = adata[col_f->source - 1];
		float flag = adata[col_rfree->source - 1];
		MaskType mask = (flag <= 0.1) ? MaskFree : MaskWork;

		long element = fft->element(h, k, l);

		fft->data[element][0] = amplitude;
		fft->data[element][1] = 0;
		fft->mask[element] = mask;

		count++;

		if (mask == MaskFree)
		{
			maskCount++;
		}
	}

	std::cout << "Loaded " << count << " reflections into"\
	" memory from " << _filename << "." << std::endl << std::endl;
	std::cout << "Counted " << maskCount << " free reflections." << std::endl;

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
