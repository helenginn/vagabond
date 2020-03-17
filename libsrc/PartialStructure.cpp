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
#include "maths.h"
#include "Diffraction.h"
#include "CSV.h"
#include "RefinementNelderMead.h"
#include "Options.h"
#include "Crystal.h"
#include "Shouter.h"

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
	
	setSolvScale(this, 0);
	setSolvBFac(this, 20);
	
	NelderMeadPtr fine;
	fine = NelderMeadPtr(new RefinementNelderMead());
	fine->setJobName("solvent_scale_fine");
	fine->setEvaluationFunction(scalePartialScore, this);
	fine->setCycles(120);
	fine->setSilent(true);
	fine->addParameter(this, getSolvScale, setSolvScale, 0.2, 1.0, "scale");
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
	VagFFTPtr fftData = _data->getFFT();
	VagFFTPtr fft = getCrystal()->getFFT();
	CSym::CCP4SPG *spg = getCrystal()->getSpaceGroup();
	mat3x3 toRecip = _partial->toRecip();
	
	vec3 nLimits = getNLimits(fft, _partial);

	std::vector<double> fData, fModel;
	double adjB = _solvBFac;
	double adjS = _solvScale;
	if (adjS < 0) { adjS = 0; }
	
	double swapped = 0;
	double lowres = 0;
	
	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				if (score)
				{
					if (!fftData->withinBounds(i, j, k))
					{
						continue;
					}

					int m, n, o;
					int asu = CSym::ccp4spg_put_in_asu(spg, i, j, k, 
					                                   &m, &n, &o);

					if (!(m == i && n == j && o == k))
					{
						continue;
					}
				}

				long fi = fftData->element(i, j, k);
				
				double test = fftData->getReal(fi);
				if (test != test && score) continue;
				
				float ref = sqrt(fftData->getIntensity(fi));
				if (score && ref != ref) continue;
				
				long nModel = fft->element(i, j, k);
				long nPart = _partial->element(i, j, k);
				float realProtein = fft->getReal(nModel);
				float imagProtein = fft->getImag(nModel);
				float realPartial = _partial->getReal(nPart);
				float imagPartial = _partial->getImag(nPart);
				
				if (realProtein != realProtein ||
				    realPartial != realPartial)
				{
					continue;
				}

				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(toRecip, &ijk);

				double length = vec3_length(ijk);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-adjB / four_d_sq);


				if (score)
				{
					bool isFree;
					isFree = fftData->getScratchComponent(fi, 0, 0) < 0.5;
					if (isFree) continue;
				}


				realPartial *= adjS * bFacMod;
				imagPartial *= adjS * bFacMod;
				
				float real = realProtein + realPartial;
				float imag = imagProtein + imagPartial;
				
				if (score)
				{
					float amp = sqrt(real * real + imag * imag);

					fData.push_back(ref);
					fModel.push_back(amp);
				}
				else
				{
					_partial->setComponent(nPart, 0, realPartial);
					_partial->setComponent(nPart, 1, imagPartial);

					if (real != real || imag != imag)
					{
						continue;
					}

					fft->setComponent(nModel, 0, real);
					fft->setComponent(nModel, 1, imag);
				}
			}
		}
	}
	
	if (!score)
	{
		return 0;
	}

	double correl = correlation(fData, fModel);
	
	/* Penalty for having a negative B factor */
	if (adjB < 0)
	{
		adjB /= 100;
		correl *= exp(adjB);
	}
	
	return -correl;
}
