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

#include "Correl.h"
#include "Crystal.h"
#include "Options.h"
#include "Diffraction.h"
#include "fftw3d.h"
#include "CSV.h"
#include <iomanip>

void Correl::setCrystalAndData(CrystalPtr crystal, DiffractionPtr data)
{
	_crystal = crystal;
	_fft = crystal->getFFT();
	_data = data->getFFT();
	_diff = data;
}

void Correl::setupShells()
{
	double minRes = Options::minRes();
	makeShells(&_shells, minRes, _crystal->getMaxResolution(_diff));
}

void Correl::localCC()
{
	setupShells();

	double xLimit = _fft->nx / 3;
	double yLimit = _fft->ny / 3;
	double zLimit = _fft->nz / 3;
	
	CSym::CCP4SPG *spg = _crystal->getSpaceGroup();
	
	const int streak = 4;
	std::vector<double> xs;
	std::vector<double> ys;
	int num = 0;
	mat3x3 r2f = _fft->getReal2Frac();
	mat3x3 tmp = mat3x3_transpose(r2f);

	for (int i = -xLimit; i < xLimit; i++)
	{
		for (int j = -yLimit; j < yLimit; j++)
		{
			for (int k = -zLimit; k < zLimit; k++)
			{
				int _i = 0; int _j = 0; int _k = 0;
				vec3 ijk = make_vec3(i, j, k);
				CSym::ccp4spg_put_in_asu(spg, i, j, k, &_i, &_j, &_k);
				
				if (!(i == _i && j == _j && k == _k))
				{
					continue;
				}

				int isFree = (_fft->getMask(_i, _j, _k) == 0);
				
				xs.clear(); ys.clear();
				
				for (int g = 0; g < 3; g++)
				{
					for (int h = -streak; h <= streak; h++)
					{
						int hi = (g % 2 == 0) * h;
						int hj = (g % 2 == 1) * h;
						int hk = (g % 2 == 2) * h;
						
						/* Only include middle value once */
						if (g > 0 && h == 0)
						{
							continue;
						}

						double int_data = _data->getReal(_i + hi, _j + hj, 
						                                      _k + hk);
						double int_model = _fft->getIntensity(i + hi, j + hj, 
						                                      k + hk);

						if (int_data != int_data || int_model != int_model)
						{
							continue;
						}

						xs.push_back(int_data);
						ys.push_back(int_model);
						
					}
				}
				
				if (xs.size() < 7)
				{
					continue;
				}

				double correl = correlation(xs, ys);

				int index = -1;

				mat3x3_mult_vec(tmp, &ijk);
				double length = vec3_length(ijk);
				double real_space = 1 / length;
				
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
				
				double sigma = _data->getImaginary(_i, _j, _k);

				if (sigma != sigma)
				{
					continue;
				}

				double iobs = _data->getReal(_i, _j, _k);
				double imodel = _fft->getIntensity(_i, _j, _k);
				
				double weight = 1 / (sigma * sigma);

				_shells[index].aveFo += correl * weight;
				_shells[index].count += weight;
				
				_shells[index].work1.push_back(sqrt(iobs));
				_shells[index].work2.push_back(sqrt(imodel));
			}
		}
	}
	
	CSVPtr csv = CSVPtr(new CSV(4, "min", "max", "local", "overall"));
	double sum = 0; double weights = 0;
	
	for (int i = 0; i < _shells.size(); i++)
	{
		_shells[i].aveFo /= _shells[i].count;
		sum += _shells[i].aveFo * _shells[i].count;
		weights += _shells[i].count;
		
		double correl = correlation(_shells[i].work1, _shells[i].work2);
		
		csv->addEntry(4, _shells[i].minRes, _shells[i].maxRes,
		              _shells[i].aveFo, correl);
	}
	
	int cycle = _crystal->getCycleNum();
	csv->setSubDirectory("correlation_plots");
	csv->writeToFile("shell_by_shell_" + i_to_str(cycle) + ".csv");
	
	double overall = sum / weights;
	std::cout << "Overall local CC: " << std::setprecision(4) << overall << std::endl;
	
	_crystal->setLastLocalCC(overall);
}
