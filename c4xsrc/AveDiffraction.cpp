// cluster4x
// Copyright (C) 2019 Helen Ginn
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

#include "AveDiffraction.h"
#include "Group.h"
#include "MtzFile.h"
#include "MtzFFT.h"
#include <libsrc/Crystal.h>
#include <libsrc/FileReader.h>

AveDiffraction::AveDiffraction(Group *group, double maxRes) : Average(group)
{
	_maxRes = maxRes;
	_origGroup = group;
}

AveDiffraction::~AveDiffraction()
{

}

void AveDiffraction::calculate()
{
	VagFFTPtr first = _mtzs[0];
	_fft = VagFFTPtr(new VagFFT(*first));
	_fft->wipe();
	std::vector<double> ucs = std::vector<double>(6, 0);
	int count = 0;

	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		MtzFFTPtr current = _mtzs[i];
		if (current->getMtzFile()->isDead())
		{
			continue;
		}
		
		std::vector<double> uc = current->getUnitCell();
		
		for (size_t j = 0; j < uc.size(); j++)
		{
			ucs[j] += uc[j];
		}
		
		count++;

		vec3 nLimits = getNLimits(_fft, current);

		for (int k = -nLimits.z; k < nLimits.z; k++)
		{
			for (int j = -nLimits.y; j < nLimits.y; j++)
			{
				for (int i = -nLimits.x; i < nLimits.x; i++)
				{
					double res = current->resolution(i, j, k);
					if (res < _maxRes)
					{
						continue;
					}

					double amp = current->getReal(i, j, k);
					long ele = _fft->element(i, j, k);
					if (amp == amp)
					{
						_fft->addToReal(ele, amp);
						_fft->addToImag(ele, 1);
					}
				}
			}
		}
	}
	
	for (size_t i = 0; i < ucs.size(); i++)
	{
		ucs[i] /= (double)count;
	}

	_fft->setUnitCell(ucs);
	
	for (long i = 0; i < _fft->nn(); i++)
	{
		double imag = _fft->getImag(i);
		double real = _fft->getReal(i);

		real /= imag;

		_fft->setReal(i, real);
	}
	
//	scaleIndividuals(_origGroup);
}

double AveDiffraction::findCorrelation(MtzFFTPtr one, MtzFFTPtr two)
{
	CorrelData cd = empty_CD();

	for (int i = 0; i < one->nn(); i++)
	{
		double ave = _fft->getReal(i);
		double mine = one->getReal(i);
		double yours = two->getReal(i);

		double x = mine - ave;
		double y = yours - ave;

		add_to_CD(&cd, x, y);
	}

	double cc = evaluate_CD(cd);
	if (cc < 0) return 0;
	return cc;
}

void AveDiffraction::findIntercorrelations(Group *other, double **svd)
{
	scaleIndividuals(other);

	Average::findIntercorrelations(other, svd);
}

void AveDiffraction::scaleIndividualMtz(MtzFFTPtr mtz)
{
	std::vector<ShellInfo> _shells;
	makeShells(&_shells, 0, _maxRes);

	vec3 nLimits = getNLimits(_fft, mtz);

	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = mtz->resolution(i, j, k);
				if (res < _maxRes)
				{
					continue;
				}

				double amp = mtz->getReal(i, j, k);
				double ref = _fft->getReal(i, j, k);

				int s = findShell(_shells, res);
				
				if (s < 0)
				{
					continue;
				}

				_shells[s].work1.push_back(amp);
				_shells[s].work2.push_back(ref);
			}
		}
	}
	
	for (size_t i = 0; i < _shells.size(); i++)
	{
		_shells[i].scale = scale_factor_by_sum(_shells[i].work1,
		                                       _shells[i].work2);
	}
	
	for (int k = -nLimits.z; k < nLimits.z; k++)
	{
		for (int j = -nLimits.y; j < nLimits.y; j++)
		{
			for (int i = -nLimits.x; i < nLimits.x; i++)
			{
				double res = mtz->resolution(i, j, k);
				int s = findShell(_shells, res);
				
				if (s < 0)
				{
					mtz->setReal(i, j, k, 0);
					continue;
				}

				double scale = _shells[s].scale;
				double amp = mtz->getReal(i, j, k);
				amp *= scale;
				mtz->setReal(i, j, k, amp);
			}
		}
	}
}

void AveDiffraction::scaleIndividuals(Group *other)
{
	std::cout << "Scaling individual data sets..." << std::flush;
	for (size_t i = 0; i < other->mtzCount(); i++)
	{
		scaleIndividualMtz(other->getMtz(i));
	}
	
	std::cout << " ... done." << std::endl;
}


std::string AveDiffraction::unitCellDesc(VagFFTPtr fft)
{
	std::vector<double> uc = fft->getUnitCell();
	std::string ucInfo;
	ucInfo += fft->getSpaceGroup()->symbol_Hall;
	ucInfo += "; a, b, c =  ";
	for (int i = 0; i < 3; i++)
	{
		ucInfo += f_to_str(uc[i], 2) + "   ";
	}
	ucInfo += " Å; a, b, g =  ";
	for (int i = 3; i < 6; i++)
	{
		ucInfo += f_to_str(uc[i], 2) + "   ";
	}
	ucInfo += "°.";

	return ucInfo;
}

