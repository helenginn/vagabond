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

#include "PartialStructure.h"
#include "CompareFFT.h"
#include <hcsrc/maths.h>
#include "Diffraction.h"
#include "CSV.h"
#include <hcsrc/RefinementNelderMead.h>
#include "Options.h"
#include "Crystal.h"
#include "Shouter.h"

PartialStructure::PartialStructure()
{
	_compare = NULL;
}

void PartialStructure::setStructure(VagFFTPtr refPart)
{
	VagFFTPtr vFFT = getCrystal()->getFFT();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();

	_partial = VagFFTPtr(new VagFFT(*vFFT));
	_partial->wipe();

	vec3 nLimits = getNLimits(refPart, _partial);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				int _h, _k, _l;
				int sym = CSym::ccp4spg_put_in_asu(spg, i, j, k, 
				                                   &_h, &_k, &_l);

				long index = _partial->element(i, j, k);
				long asuIndex = _partial->element(_h, _k, _l);
				long pIndex = refPart->element(_h, _k, _l);

				double x = refPart->getReal(pIndex);
				double y = refPart->getImag(pIndex);

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

				_partial->setComponent(index, 0, x);
				_partial->setComponent(index, 1, y);
			}
		}
	}
}

void PartialStructure::reportScale()
{
	if (!getCrystal()->isSilent())
	{
		std::cout << "   Partial scale: " << getSolvScale(this) << ", "
		"B factor: " << getSolvBFac(this) << std::endl;
	}
}

void PartialStructure::scalePartialStructure()
{
	if (!_data)
	{
		shout_at_helen("Need diffraction data to scale solvent");
	}

	std::string cycle = i_to_str(getCrystal()->getCycleNum());

	VagFFTPtr fft = getCrystal()->getFFT();
	VagFFTPtr fftData = _data->getFFT();

	if (_compare != NULL)
	{
		delete _compare;
		_compare = NULL;
	}

	_compare = new CompareFFT();
	_compare->setPrimaryFFT(fftData);
	_compare->setSecondaryFFT(fft);
	_compare->setTertiaryFFT(_partial);
	_compare->setupResolutions(true);
//	_compare->setResolutionCutoff(2.0);
	_compare->prepare();

	setSolvScale(this, 0);
	setSolvBFac(this, 40);

	NelderMeadPtr fine;
	fine = NelderMeadPtr(new RefinementNelderMead());
	fine->setJobName("solvent_scale_fine");
	fine->setEvaluationFunction(scalePartialScore, this);
	fine->setCycles(120);
	fine->setSilent(true);
	fine->addParameter(this, getSolvScale, setSolvScale, 0.2, 0.0001, "scale");
	fine->addParameter(this, getSolvBFac, setSolvBFac, 40, 0.1, "bfac");
	fine->refine();

	if (_solvScale < 0)
	{
		_solvScale = 0;
	}

	/** Now add this into the FFT */

	scaleOrScore(false);
}

double PartialStructure::scalePartialScore(void *object)
{
	return static_cast<PartialStructure *>(object)->scaleAndAddPartialScore();
}

double PartialStructure::scaleAndAddPartialScore()
{
	return scaleOrScore(true);

}

double PartialStructure::scaleOrScore(bool score)
{
	VagFFTPtr observed = _data->getFFT();
	VagFFTPtr calculated = getCrystal()->getFFT();

	CorrelData cd = empty_CD();

	double adjB = _solvBFac;
	double adjS = _solvScale;
	if (adjS < 0) { adjS = 0; }
	
	bool found = false;
	
	long num = score ? _compare->pairCount() : _compare->allPairCount();

	for (size_t i = 0; i < num; i++)
	{
		CompareFFT::FFTPair &pair = (score ? _compare->pair(i) : 
		                             _compare->allPair(i));


		long fi = pair.id1;
		float ref = sqrt(pair.data1[0]* pair.data1[0] + 
		                 pair.data1[1] * pair.data1[1]);

		long nModel = pair.id2;
		float realProtein = calculated->getReal(nModel);
		float imagProtein = calculated->getImag(nModel);

		float realPartial = pair.data3[0];
		float imagPartial = pair.data3[1];

		double d = pair.resolution;
		double four_d_sq = (4 * d * d);
		double bFacMod = exp(-adjB / four_d_sq);

		realPartial *= adjS * bFacMod;
		imagPartial *= adjS * bFacMod;
		
		float real = realProtein + realPartial;
		float imag = imagProtein + imagPartial;

		if (score)
		{
			float amp = sqrt(real * real + imag * imag);
			if (real != real || imag != imag)
			{
				continue;
			}
			
			if (i > 4000 && !found)
			{
				found = true;
			}

			add_to_CD(&cd, ref, amp);
		}
		else
		{
			long nPart = pair.id3;
			_partial->setComponent(nPart, 0, realPartial);
			_partial->setComponent(nPart, 1, imagPartial);

			if (real != real || imag != imag)
			{
				continue;
			}

			calculated->setComponent(nModel, 0, real);
			calculated->setComponent(nModel, 1, imag);
		}
	}
	
	if (!score)
	{
		return 0;
	}

	double correl = evaluate_CD(cd);

	/* Penalty for having a negative B factor */
	if (adjB < 0)
	{
		adjB /= 100;
		correl *= exp(adjB);
	}

	return -correl;
}
