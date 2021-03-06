// Vagabond
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

#include "FFT.h"
#include "Options.h"
#include "CSV.h"
#include "Element.h"
#include "Atom.h"
#include "Timer.h"
#include "Anisotropicator.h"
#include "ExplicitModel.h"
#include <cstring>
#include <float.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#include "libccp4/cmtzlib.h"
#include "libccp4/csymlib.h"
#include "libccp4/ccp4_spg.h"
#include "libccp4/ccp4_general.h"

std::vector<FFTDim *> VagFFT::_dimensions;

void checkDataPtr(fftwf_complex *data, long nn)
{
	if (!data)
	{
		std::cout << "Error: likely exceeded available memory." << std::endl;
		std::cout << "Malloc failed for VagFFT, nn = " << nn << std::endl;
		exit(1);
	}
}


VagFFT::VagFFT(VagFFT &fft, int scratch)
{
	_nx = fft._nx;
	_ny = fft._ny;
	_nz = fft._nz;
	_nele = fft._nele;
	_nscratch = fft._nscratch;
	_activeScratch = fft._activeScratch;

	_bFactor = fft._bFactor;
	_nn = fft._nn;
	_status = fft._status;
	_stride = fft._stride;
	_total = fft._total;
	_origin = fft._origin;
	_toReal = fft._toReal;
	_toRecip = fft._toRecip;
	_recipTrans = fft._recipTrans;
	_recipBasis = fft._recipBasis;
	_realBasis = fft._realBasis;
	_setMatrices = fft._setMatrices;
	_elements = fft._elements;
	_elementMap = fft._elementMap;
	_spg = fft._spg;
	_lowResMode = fft._lowResMode;
	_unitCell = fft._unitCell;
	_data = NULL;
	
	if (scratch >= 0 && scratch != fft._nscratch)
	{
		_nscratch = scratch; 
		_stride = 2 * (_nele) + 1 + _nscratch;
		_total = _nn * _stride * sizeof(fftwf_complex);
		_data = (fftwf_complex *)fftwf_malloc(_total);
		int start = finalIndex(0);
		_lastData = &_data[start];
		checkDataPtr(_data, _nn);
		memset(_data, '\0', _total);
		_myDims = 0;
	}
	else
	{
		_myDims = fft._myDims;
		_data = (fftwf_complex *)fftwf_malloc(_total);
		checkDataPtr(_data, _nn);
		int start = finalIndex(0);
		_lastData = &_data[start];
		memcpy(_data, fft._data, _total);
	}
}

void VagFFT::prepareShortScratch()
{
	_shortScratch = (float *)fftwf_malloc(_nn * sizeof(float));
	memset(_shortScratch, '\0', _nn * sizeof(float));
}

VagFFT::~VagFFT()
{
	if (_data != NULL)
	{
		fftwf_free(_data);
	}
	
	_data = NULL;
}

VagFFT::VagFFT(int nx, int ny, int nz, int nele, int scratches)
{
	_bFactor = 0;
	_shortScratch = NULL;
	_nx = nx;
	_ny = ny;
	_nz = nz;
	
	if (_nx % 2 == 1) _nx -= 1;
	if (_ny % 2 == 1) _ny -= 1;
	if (_nz % 2 == 1) _nz -= 1;

	_myDims = NULL;
	_spg = CSym::ccp4spg_load_by_ccp4_num(1);
	_lowResMode = Options::getLowResMode();

	_nele = nele;
	_nscratch = scratches;
	_activeScratch = 0;
	_nn = _nx * _ny * _nz;
	_status = FFTEmpty;

	_stride = 2 * (_nele) + 1 + _nscratch;
	_total = _nn * _stride * sizeof(fftwf_complex);
	_data = (fftwf_complex *)fftwf_malloc(_total);
	checkDataPtr(_data, _nn);

	int start = finalIndex(0);
	_lastData = &_data[start];

	memset(_data, 0, _total);

	_origin = empty_vec3();
	_toReal = make_mat3x3();
	_toRecip = make_mat3x3();
	_recipTrans = make_mat3x3();
	_realBasis = make_mat3x3();
	_recipBasis = make_mat3x3();
	_setMatrices = false;
}

void VagFFT::addElement(ElementPtr ele)
{
	_elementMap[ele] = _elements.size();
	_elements.push_back(ele);
}

void VagFFT::printStatus()
{
	std::cout << "Current status: " << std::flush;
	
	switch (_status)
	{
		case FFTEmpty:
		std::cout << "FFTEmpty";
		break;
		case FFTRealSpace:
		std::cout << "FFTRealSpace";
		break;
		case FFTReciprocalSpace:
		std::cout << "FFTReciprocalSpace";
		break;
		case FFTSeparateAtoms:
		std::cout << "FFTSeperateAtoms";
		break;
		default:
		std::cout << "Unknown";
		break;
	}
	
	std::cout << std::endl;
	
	std::cout << "Dimensions: " << _nx << " " << _ny <<
	" " << _nz << std::endl;
	
	for (int i = -1; i < _nscratch; i++)
	{
		std::cout << "Amplitude for scratch " << i <<
		" " << std::setprecision(6) << sumAmp(i) << std::endl;
	}
	
	std::cout << std::endl;
}

bool VagFFT::sanityCheck()
{
	if (_nx <= 0 || _ny <= 0 || _nz <= 0)
	{
		std::cout << "Illogical VagFFT dimensions: " << _nx 
		<< " " << _ny << " " << _nz << std::endl;
		return false;
	}
	
	if (_nele != _elements.size())
	{
		std::cout << "Element number and assigned element count do"\
		" not match, " << _nele << " vs. " << _elements.size() << std::endl;
		return false;
	}
	
	if (!_setMatrices)
	{
		std::cout << "Need to set up real / reciprocal change of bases" 
		<< std::endl;
		return false;
	}
	
	return true;
}

void VagFFT::wipe()
{
	memset(_data, 0, _total);
	_status = FFTEmpty;
}

int scoreN(int n)
{
	int total = 0;
	for (int i = 3; i < n / 2; i++)
	{
		total += (n % i == 0 ? 1 : 0);
	}
	
	return total;
}

void adjustN(int *n)
{
	int nMin = *n;
	int nMax = *n * 1.1;
	int best = scoreN(*n);

	for (int i = nMin; i < nMax; i++)
	{
		if (i % 2 == 1)
		{
			continue;
		}

		int score = scoreN(i);
		if (score > best)
		{
			best = score;
			*n = i;
		}
	}
}

void VagFFT::adjustNs()
{

	double xtmp = _nx;
	double ytmp = _ny;
	double ztmp = _nz;
	adjustN(&_nx);
	adjustN(&_ny);
	adjustN(&_nz);
	
	_nn = _nx * _ny * _nz;
	
	_total = _nn * _stride * sizeof(fftwf_complex);
	
	if (_data != NULL)
	{
		fftwf_free(_data);
		_data = (fftwf_complex *)fftwf_malloc(_total);
		checkDataPtr(_data, _nn);
		memset(_data, 0, _total);
		int start = finalIndex(0);
		_lastData = &_data[start];
	}

	_realBasis.vals[0] *= _nx / xtmp;
	_realBasis.vals[4] *= _ny / ytmp;
	_realBasis.vals[8] *= _nz / ztmp;
}

void VagFFT::makePlans()
{
	if (!sanityCheck())
	{
		std::cout << "Sanity check failed" << std::endl;
		return;
	}

	/* do we already have this planned dimension? */
	for (int i = 0; i < _dimensions.size(); i++)
	{
		if (_dimensions[i]->nx == _nx
		    && _dimensions[i]->ny == _ny
		    && _dimensions[i]->nz == _nz
		    && _dimensions[i]->nele == _nele
		    && _dimensions[i]->nscratch == _nscratch)
		{
			_myDims = _dimensions[i];
			return;
		}
	}
	
	Timer timer;

	FFTDim *dims = (FFTDim *)malloc(sizeof(FFTDim));
	dims->nx = _nx; dims->ny = _ny; dims->nz = _nz; dims->nele = _nele;
	dims->nscratch = _nscratch;
	_dimensions.push_back(dims);
	_myDims = _dimensions.back();
	
	int ns[3] = {_nz, _ny, _nx};

	/* destroying input should do nothing though */
	unsigned many_flags = FFTW_MEASURE;

	fftwf_plan plan = fftwf_plan_many_dft(3, ns, _nele, _data, NULL, 
	                                      _stride, 2, _data, 
	                                      NULL, _stride, 2, 
	                                      1, many_flags); 
	_myDims->atom_real_to_recip = plan;

	plan = fftwf_plan_many_dft(3, ns, _nele, _data, NULL, _stride, 
	                                      2, _data, 
	                                      NULL, _stride, 2, 
	                                      -1, many_flags); 
	_myDims->atom_recip_to_real = plan;

	plan = fftwf_plan_many_dft(3, ns, 1, _lastData, NULL, /*null=embed*/
	                           _stride, 2, _lastData, NULL, 
	                           _stride, 2, 1, many_flags); 

	_myDims->real_to_recip = plan;

	plan = fftwf_plan_many_dft(3, ns, 1, _lastData, NULL, 
	                           _stride, 2, _lastData, NULL,
	                           _stride, 2, -1, many_flags); 

	_myDims->recip_to_real = plan;
}

void VagFFT::multiplyDotty(float val)
{
	for (int i = 0; i < _nn; i++)
	{
		for (int j = 0; j < _nele; j++)
		{
			long dotty_index = dottyIndex(i, j);
			_data[dotty_index][0] *= val;
			_data[dotty_index][1] *= val;
		}
	}
}

void VagFFT::multiplyFinal(float val)
{
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][0] *= val;
		_data[final_index][1] *= val;
	}
}


void VagFFT::setupElements(bool wipe)
{
	if (_elements.size() == 0 && _nele == 0)
	{
		return;
	}
	
	for (int z = -_nz / 2; z <= _nz / 2; z++)
	{
		for (int y = -_ny / 2; y <= _ny / 2; y++)
		{
			for (int x = -_nx / 2; x <= _nx / 2; x++)
			{
				long ele = element(x, y, z);
				
				vec3 real = make_vec3(x, y, z);
				mat3x3_mult_vec(_toRecip, &real);

				double length = vec3_length(real);
				double d = 1 / length;
				double four_d_sq = (4 * d * d);
				double bFacMod = exp(-_bFactor / four_d_sq);

				vec3_mult(&real, 0.5);

				for (int j = 0; j < _nele; j++)
				{
					long ele_index = eleIndex(ele, j);
					
					ElementPtr elem = _elements[j];
					double val = Element::getVoxelValue(&*elem, real.x, 
					                                    real.y, real.z);
					val *= bFacMod;

					_data[ele_index][0] = val;
					_data[ele_index][1] = 0;
					
					if (wipe)
					{
						long dotty_index = ele_index - 1;
						_data[dotty_index][0] = 0; 
						_data[dotty_index][1] = 0;
					}
				}
			}
		}
	}
}

void VagFFT::copyScratchElementToPosition(ElementPtr ele)
{
	int column = whichColumn(ele);

	for (int i = 0; i < _nn; i++)
	{
		_data[i * _stride + column][0] = _shortScratch[i]; 
	}

	memset(_shortScratch, '\0', _nn * sizeof(float));
}

void VagFFT::prepareAtomSpace()
{
	if (_status == FFTEmpty)
	{
		setupElements();
		prepareShortScratch();
	}

	_status = FFTSeparateAtoms;
}

void VagFFT::separateAtomTransform()
{
	fftwf_execute_dft(_myDims->atom_real_to_recip, _data, _data);
	
	/* go through and multiply transform of atoms by element factor */
	for (int i = 0; i < _nn; i++)
	{
		long final_idx = finalIndex(i);
		_data[final_idx][0] = 0;
		_data[final_idx][1] = 0;

		for (int j = 0; j < _nele; j++)
		{
			long dotty_idx = dottyIndex(i, j);
			long ele_idx = dotty_idx + 1;

			_data[final_idx][0] += _data[dotty_idx][0] * _data[ele_idx][0];
			_data[final_idx][1] += _data[dotty_idx][1] * _data[ele_idx][0];
			
			_data[dotty_idx][0] = 0;
			_data[dotty_idx][1] = 0;
		}
	}

	_status = FFTReciprocalSpace;
}

void VagFFT::fft(FFTTransform transform)
{
	if (_myDims == NULL)
	{
		std::cout << "Attempting to transform FFT but not set up any"
		" plans for doing so." << std::endl;
		exit(1);
	}

	if (transform == FFTAtomsToReciprocal || transform == FFTAtomsToReal)
	{
		if (!(_status == FFTSeparateAtoms || _status == FFTEmpty))
		{
			std::cout << "Trying to transform individual atoms when status is "
			<< "not FFTSeparateAtoms" << std::endl;
			printStatus();
			exit(1);
		}

		separateAtomTransform();

		if (transform == FFTAtomsToReal)
		{
			fftwf_execute_dft(_myDims->recip_to_real, _lastData, _lastData);
			multiplyFinal(1 / (double)_nn);
			_status = FFTRealSpace;
		}
	}

	if (transform == FFTRealToReciprocal)
	{
		if (!(_status == FFTRealSpace || _status == FFTEmpty))
		{
			std::cout << "Asking for transform Real to Reciprocal, but "\
			"FFT status is not Real Space" << std::endl;
			printStatus();
			exit(1);
		}

		fftwf_execute_dft(_myDims->real_to_recip, _lastData, _lastData);
		_status = FFTReciprocalSpace;
	}

	if (transform == FFTReciprocalToReal)
	{
		if (!(_status == FFTReciprocalSpace || _status == FFTEmpty))
		{
			std::cout << "Asking for transform Reciprocal to Real, but "\
			"FFT status is not Reciprocal Space" << std::endl;
			printStatus();
			exit(1);
		}

		fftwf_execute_dft(_myDims->recip_to_real, _lastData, _lastData);
		multiplyFinal(1 / (double)_nn);
		_status = FFTRealSpace;
	}

}

void VagFFT::setUnitCell(std::vector<double> dims)
{
	_toReal = mat3x3_from_unit_cell(dims[0], dims[1], dims[2], dims[3],
	                                dims[4], dims[5]);
	_toRecip = mat3x3_inverse(_toReal);
	_recipTrans = mat3x3_transpose(_toRecip);

	_realBasis = _toReal;
	mat3x3_scale(&_realBasis, 1 / (double)_nx, 1 / (double)_ny, 
	             1 / (double)_nz);
	_recipBasis = mat3x3_inverse(_realBasis);
	_setMatrices = true;
	_unitCell = dims;
}

void VagFFT::addImplicitAtom(AtomPtr atom)
{
	mat3x3 tensor = atom->getModel()->getRealSpaceTensor();

	Anisotropicator aniso;
	aniso.setTensor(tensor);
	mat3x3 ellipsoid = aniso.basis();
	
	/* make these variances into standard deviations */
	for (int i = 0; i < 3; i++)
	{
		vec3 axis = mat3x3_axis(ellipsoid, i);
		double length = vec3_length(axis);
		
		for (int j = 0; j < 3; j++)
		{
			ellipsoid.vals[j * 3 + i] /= sqrt(length);
		}
	}
	
	vec3 centre = atom->getAbsolutePosition();
	
	const int scale = 3.0; /* no. standard devs we care about */

	vec3 maxVals = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			*(&maxVals.x + i) = std::max(*(&maxVals.x + i),
			                             fabs(ellipsoid.vals[i * 3 + j]));
		}
	}
	
	vec3_mult(&maxVals, scale);
	maxVals.x += 1; maxVals.y += 1; maxVals.z += 1;
	/* into voxel coordinates */
	mat3x3_mult_vec(_recipBasis, &maxVals);
	
	vec3_subtract_from_vec3(&centre, _origin);
	mat3x3_mult_vec(_recipBasis, &centre);
	ElementPtr ele = atom->getElement();
	int column = whichColumn(ele);
	
	double occ = atom->getModel()->getEffectiveOccupancy();
	populateImplicit(centre, maxVals, tensor, ellipsoid, occ);
}

double VagFFT::populateImplicit(vec3 centre, vec3 maxVals,
                                mat3x3 tensor, mat3x3 ellipsoid, 
                                double scale)
{
	double total = 0;
	double count = 0;
	
	const double factor = 1 / sqrt(2 * M_PI);
	double squash = factor * factor * factor;
	/* make these variances into standard deviations */
	for (int i = 0; i < 3; i++)
	{
		vec3 axis = mat3x3_axis(ellipsoid, i);
		squash /= vec3_length(axis);
	}
	
	double vol = mat3x3_volume(_realBasis);
	mat3x3 inv = mat3x3_inverse(tensor);
	
	vec3 whole = make_vec3((int)floor(centre.x), 
                           (int)floor(centre.y),
                           (int)floor(centre.z));

	vec3 rest = make_vec3(centre.x - whole.x, 
	                      centre.y - whole.y,
	                      centre.z - whole.z);

	for (int z = -maxVals.z; z <= maxVals.z + 0.5; z++)
	{
		for (int y = -maxVals.y; y <= maxVals.y + 0.5; y++)
		{
			for (int x = -maxVals.x; x <= maxVals.x + 0.5; x++)
			{
				vec3 pos = make_vec3(x, y, z);
				pos = make_vec3(x + centre.x, y + centre.y, z + centre.z);
				/* in voxel system */
				vec3 offset = vec3_subtract_vec3(pos, centre);
				/* and then to Angstroms */
				mat3x3_mult_vec(_realBasis, &offset);

				vec3 orig = offset;

				/* to aniso-U-sensitive coordinates */
				mat3x3_mult_vec(inv, &offset);

				offset.x *= orig.x;
				offset.y *= orig.y;
				offset.z *= orig.z;
				double mult = offset.x + offset.y + offset.z;
				
				double dens = exp(-mult / 2.0);
				dens *= squash * scale * vol;

				addInterpolatedToReal(x + centre.x, y + centre.y, 
				                      z + centre.z, dens);
			}
		}
	}

	return total;
}

int VagFFT::whichColumn(ElementPtr ele)
{
	if (_elementMap.count(ele))
	{
		return _elementMap[ele] * 2;
	}

	return -1;
}

void VagFFT::calculateAdjustmentVolumes()
{
	double total = 0;
	std::vector<double> ord;
	ord.push_back(_realBasis.vals[0]);
	ord.push_back(_realBasis.vals[4]);
	ord.push_back(_realBasis.vals[8]);
	std::sort(ord.begin(), ord.end(), std::greater<double>());

	_vols[0] = _realBasis.vals[0] / ord[2];
	_vols[1] = _realBasis.vals[4] / ord[2];
	_vols[2] = _realBasis.vals[8] / ord[2];
}

void VagFFT::addInterpolatedToReal(double sx, double sy, 
                                   double sz, double val)
{
	collapse(&sx, &sy, &sz);

	double ss[3] = {sx, sy, sz};

	double ls[3], rs[3];
	
	for (int i = 0; i < 3; i++)
	{
		ls[i] = (int)floor(ss[i]);
		rs[i] = ss[i] - ls[i];
		
		if (rs[i] < 0.5)
		{
			rs[i] *= _vols[i];
		}
		else
		{
			rs[i] = 1 - (1 - rs[i]) * _vols[i];
		}
	}
	
	sx = ls[0] + rs[0];
	sy = ls[1] + rs[1];
	sz = ls[2] + rs[2];
	
	double xProps[2];
	double yProps[2];
	double zProps[2];

	xProps[0] = std::max(1 + ls[0] - sx, 0.);
	xProps[1] = std::max(sx - ls[0], 0.);

	yProps[0] = std::max(1 + ls[1] - sy, 0.);
	yProps[1] = std::max(sy - ls[1], 0.);

	zProps[0] = std::max(1 + ls[2] - sz, 0.);
	zProps[1] = std::max(sz - ls[2], 0.);
	
	for (int r = 0; r < 2; r++)
	{
		for (int q = 0; q < 2; q++)
		{
			for (int p = 0; p < 2; p++)
			{
				int sx1 = ls[0] + p;
				int sy1 = ls[1] + q;
				int sz1 = ls[2] + r;

				long index = element(sx1, sy1, sz1);
				double prop = xProps[p] * yProps[q] * zProps[r];
				
				_shortScratch[index] += prop * val;
			}	
		}
	}
}

double VagFFT::getPhase(long i)
{
	long index = finalIndex(i);

	double degrees = atan2(_data[index][1], _data[index][0]) * 180 / M_PI;

	while (degrees >= 360) degrees -= 360;

	while (degrees < 0) degrees += 360;

	return degrees;

}

double VagFFT::getPhase(int x, int y, int z)
{
	long pre_index = element(x, y, z);
	long index = finalIndex(pre_index);

	double degrees = atan2(_data[index][1], _data[index][0]) * 180 / M_PI;

	while (degrees >= 360) degrees -= 360;

	while (degrees < 0) degrees += 360;

	return degrees;

}

double VagFFT::getAmplitude(ElementPtr ele, int x, int y, int z)
{
	int col = whichColumn(ele);
	long index = element(x, y, z) * _stride + col + 1;

	double val =  (_data[index][0] * _data[index][0] + 
	               _data[index][1] * _data[index][1]);
	
	return sqrt(val);
}

void VagFFT::drawSlice(int zVal, std::string filename)
{
	CSVPtr csv = CSVPtr(new CSV(4, "x", "y", "z", "val"));
	csv->setSubDirectory("slices");

	for (int j = 0; j < ny(); j++)
	{
		for (int i = 0; i < nx(); i++)
		{
			if (_status == FFTReciprocalSpace)
			{
				if (i == 0 && j == 0)
				{
					continue;
				}
			}
			int index = element(i, j, zVal);
			double s = getAmplitude(index);
			
			csv->addEntry(4, (double)i, (double)j, (double)zVal, s);
		}
	}

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = filename;
	plotMap["height"] = "800";
	plotMap["width"] = "800";
	plotMap["xHeader0"] = "x";
	plotMap["yHeader0"] = "y";
	plotMap["zHeader0"] = "val";

	plotMap["xTitle0"] = "a dim";
	plotMap["yTitle0"] = "b dim";
	plotMap["style0"] = "heatmap";
	plotMap["stride"] = i_to_str(nx());

	csv->writeToFile(filename + ".csv");
	csv->plotPNG(plotMap);
}

void VagFFT::printSlice(double zVal, double scale)
{
	sanityCheck();

	if (scale < 0)
	{
		scale = averageAll() * 5;
	}
	
	zVal *= _nz;
	
	for (int j = 0; j < _ny; j++)
	{
		std::cout << "| ";
		for (int i = 0; i < _nx; i++)
		{
			std::string symbol = " ";
			double value = getAmplitude(i, j, zVal);

			if (value > 0.01 * scale) symbol = ".";
			if (value > 0.02 * scale) symbol = ":";
			if (value > 0.04 * scale) symbol = "\"";
			if (value > 0.08 * scale) symbol = "*";
			if (value > 0.16 * scale) symbol = "x";
			if (value > 0.32 * scale) symbol = "H";
			if (value > 0.64 * scale) symbol = "#";
			if (value > 1.00 * scale) symbol = "@";

			std::cout << symbol;
		}

		std::cout << " |" << std::endl;
	}
	std::cout << std::endl;

}

void VagFFT::addExplicitAtom(AtomPtr atom)
{
	std::vector<BondSample> positions;
	positions = atom->getExplicitModel()->getFinalPositions();

	for (int i = 0; i < positions.size(); i++)
	{
		vec3_subtract_from_vec3(&positions[i].start, _origin);
		/* orthogonal simplification */
		mat3x3_mult_vec(_recipBasis, &positions[i].start);

		double occ = positions[i].occupancy;

		addInterpolatedToReal(positions[i].start.x, 
		                      positions[i].start.y, 
		                      positions[i].start.z, occ);
	}
}

void VagFFT::addAtom(AtomPtr atom)
{
	if (atom->getModel()->hasExplicitPositions())
	{
		addExplicitAtom(atom);
	}
	else
	{
		addImplicitAtom(atom);
	}
}

/* 11-point interpolation - attempted transcription from Dave's function
 * from GAP */
double VagFFT::cubic_interpolate(vec3 vox000, size_t im)
{
	/* vox000 has integer and real components */
	collapse(&vox000.x, &vox000.y, &vox000.z);
	
	/* Pick out just the real components - this is faster
	 * than modf */
	vec3 uvw = make_vec3(vox000.x - floor(vox000.x),
	                     vox000.y - floor(vox000.y),
	                     vox000.z - floor(vox000.z));

	/* Extra refers to the additional index to be fished
	 * for 11-point interpolation. We already get 0 and 1. */
	int extra[3] = {-1, -1, -1};
	int central[3] = {0, 0, 0};
	int next[3] = {1, 1, 1};
	
	/* If uvw components are greater than 0.5, then flip them 
	 * make the extra index one ahead and reverse the order */
	for (int i = 0; i < 3; i++)
	{
		if (*(&uvw.x + i) > 0.5)
		{
			extra[i] = 2;
			central[i] = 1;
			next[i] = 0;
			*(&uvw.x + i) = 1 - *(&uvw.x + i);
		}
	}

	int vox000x = vox000.x + central[0];
	int vox000y = vox000.y + central[1];
	int vox000z = vox000.z + central[2];
	int vox000xm = vox000.x + next[0];
	int vox000ym = vox000.y + next[1];
	int vox000zm = vox000.z + next[2];
	int vox000xn = vox000.x + extra[0];
	int vox000yn = vox000.y + extra[1];
	int vox000zn = vox000.z + extra[2];

	collapse(&vox000x, &vox000y, &vox000z);
	collapse(&vox000xm, &vox000ym, &vox000zm);
	collapse(&vox000xn, &vox000yn, &vox000zn);

	vox000y  *= _nx;
	vox000ym *= _nx;
	vox000yn *= _nx;
	vox000z  *= _nx * _ny;
	vox000zm *= _nx * _ny;
	vox000zn *= _nx * _ny;

	long idx000 = vox000x + vox000y + vox000z;
	long idx100 = vox000xm + vox000y + vox000z;
	long idx010 = vox000x + vox000ym + vox000z;
	long idx110 = vox000xm + vox000ym + vox000z;
	long idx001 = vox000x + vox000y + vox000zm;
	long idx101 = vox000xm + vox000y + vox000zm;
	long idx011 = vox000x + vox000ym + vox000zm;
	long idx111 = vox000xm + vox000ym + vox000zm;
	
	long idxn00 = vox000xn + vox000y + vox000z;
	long idx0n0 = vox000x + vox000yn + vox000z;
	long idx00n = vox000x + vox000y + vox000zn;
	
	double u = uvw.x;
	double v = uvw.y;
	double w = uvw.z;
	
	double p000 = _data[finalIndex(idx000)][im];
	double p001 = _data[finalIndex(idx001)][im];
	double p010 = _data[finalIndex(idx010)][im];
	double p011 = _data[finalIndex(idx011)][im];
	double p100 = _data[finalIndex(idx100)][im];
	double p101 = _data[finalIndex(idx101)][im];
	double p110 = _data[finalIndex(idx110)][im];
	double p111 = _data[finalIndex(idx111)][im];
	
	double a = p100 - p000;
	double b = p010 - p000;
	double c = p110 - p010;
	double d = p101 - p001;
	
	double pn00 = _data[finalIndex(idxn00)][im];
	double p0n0 = _data[finalIndex(idx0n0)][im];
	double p00n = _data[finalIndex(idx00n)][im];

	double p8value = p000+u*(a+w*(-a+d)+v*((c-a)+w*( a-c-d-p011+p111)))
	+ v*(b+w*(-p001+p011-b))+w*(-p000+p001);
	
	double mod = (p000 - 0.5 * p100 - 0.5 * pn00) * (u - u * u);
	mod += (p000 - 0.5 * p010 - 0.5 * p0n0) * (v - v * v);
	mod += (p000 - 0.5 * p001 - 0.5 * p00n) * (w - w * w);
	
	double p11value = p8value + 0.4 * mod;

	return p11value;
}

/*  For multiplying point-wise
 *  No assumption that interpolation is not needed.
 */
double VagFFT::operation(VagFFTPtr fftCrystal, VagFFTPtr fftAtom,
                      MapScoreType mapScoreType, std::vector<CoordVal> *vals,
                      bool sameScale, bool fcOnly, int step)
{
	/* I rarely comment something so heavily but I will get confused if
	 * I don't, this time, as I can't soak the protocol into the variable
	 * names. Bear in mind the three coordinate systems:
	 * (a) Angstroms
	 * (b) Crystal voxels
	 * (c) Atom voxels */

	/* find offset from origins of two maps */
	vec3 add = vec3_subtract_vec3(fftAtom->_origin, fftCrystal->_origin);

	/* to reciprocal space coordinates */
	mat3x3_mult_vec(fftCrystal->_toRecip, &add);

	/* Bring the fractional coordinate of the atom into range 0 < frac <= 1 */
	VagFFT::collapseFrac(&add.x, &add.y, &add.z);

	/* Multiply by the relative dimensions of the crystal */
	double multX = add.x * fftCrystal->_nx;
	double multY = add.y * fftCrystal->_ny;
	double multZ = add.z * fftCrystal->_nz;

	/* Get the remainder after subtracting a number of whole voxels */
	vec3 atomOffset; // in crystal coordinates, converted later.
	atomOffset.x = fmod(multX, 1.);
	atomOffset.y = fmod(multY, 1.);
	atomOffset.z = fmod(multZ, 1.);

	/* Store the non-remainder whole voxel values for way later. */
	vec3 cornerCrystal = make_vec3((int)multX, (int)multY, (int)multZ);

	/* Prepare a matrix to convert crystal voxels into atomic voxels */
	mat3x3 crystal2AtomVox = mat3x3_mult_mat3x3(fftAtom->getRecipBasis(),
	                                            fftCrystal->getRealBasis());

	/* Prepare a matrix to convert atomic voxels into crystal voxels */
	mat3x3 atomVox2Crystal = mat3x3_mult_mat3x3(fftCrystal->getRecipBasis(),
	                                            fftAtom->getRealBasis());

	/* Apply this offset and reverse it. 
	 * This small offset must be added to all future atomic coordinates prior
	 * to interpolation. This will be subtracted from an atom voxel to find
	 * the relevant crystal voxel, and should therefore be in atom voxels */
	mat3x3_mult_vec(crystal2AtomVox, &atomOffset);
	vec3_mult(&atomOffset, -1);

	/* Determine bounding box - 9th Dec 2017 */
	vec3 minAtom = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 maxAtom = make_vec3(0, 0, 0);

	for (int k = 0; k <= fftAtom->_nz; k += fftAtom->_nz)
	{
		for (int j = 0; j <= fftAtom->_ny; j += fftAtom->_ny)
		{
			for (int i = 0; i <= fftAtom->_nx; i += fftAtom->_nx)
			{
				vec3 atomCorner = make_vec3(i, j, k);
				vec3 toCrystal = mat3x3_mult_vec(atomVox2Crystal, atomCorner);
				toCrystal.x = (int)toCrystal.x;
				toCrystal.y = (int)toCrystal.y;
				toCrystal.z = (int)toCrystal.z;

				vec3_min_each(&minAtom, toCrystal);
				vec3_max_each(&maxAtom, toCrystal);
			}
		}
	}
	
	long count = 0;

	/* min/maxAtoms are in crystal coordinates.*/
	for (double k = minAtom.z; k < maxAtom.z; k += 1)
	{
		for (double j = minAtom.y; j < maxAtom.y; j += 1)
		{
			for (double i = minAtom.x; i < maxAtom.x; i += 1)
			{
				if (int(i + j + k) % step > 0)
				{
					continue;
				}

				/* Position currently in crystal coords - change to atom. */
				vec3 crystalPos = make_vec3(i, j, k);
				vec3 atomPos = crystalPos;

				if (!sameScale)
				{
					/* now we convert this into atomic voxels */
					mat3x3_mult_vec(crystal2AtomVox, &atomPos);

					if (atomPos.x < 0 || atomPos.y < 0 || atomPos.z < 0)
					{
						continue;
					}

					if (atomPos.x > fftAtom->_nx || atomPos.y > fftAtom->_ny
					    || atomPos.z > fftAtom->_nz)
					{
						continue;
					}
				}

				/* Now we must find the relative crystal voxel to write this
				 * density value to, given that the atom was wrapped around
				 * the origin (center). This should work regardless of odd/
				 * even dimension lengths. */

				/* We add the tiny offset which resulted from the atom
				 * falling between two voxels, in crystal voxels */
				vec3_add_to_vec3(&atomPos, atomOffset);
				
				/* Find the interpolated value which atomPos falls on */
				double atomReal = 0; double atomImag = 0;
				fftAtom->collapse(&atomPos.x, &atomPos.y, &atomPos.z);
				atomReal = fftAtom->cubic_interpolate(atomPos, 0);

				if (vals != NULL)
				{
					atomImag = fftAtom->cubic_interpolate(atomPos, 1);
				}

				/* We add the crystal offset so we don't end up with thousands
				 * of atoms at the very centre of our map */
				vec3 cVox = empty_vec3();
				long cIndex = 0;
				
				if (!fcOnly)
				{
					/* this is where we add the offset from the corner */
					cVox = vec3_add_vec3(crystalPos, cornerCrystal);

					fftCrystal->collapse(&cVox.x, &cVox.y, &cVox.z);

					/* Get the index of this final crystal voxel.
 * 					+0.5 ensures no rounding errors with floating points */
					cIndex = fftCrystal->element(cVox.x + 0.5,
					                             cVox.y + 0.5,
					                             cVox.z + 0.5);
				}

				if (mapScoreType == MapScoreTypeCorrel)
				{
					if (!fcOnly)
					{
						// We do NOT need to interpolate //
						double realCryst = fftCrystal->getReal(cIndex);

						if (vals)
						{
							CoordVal val;
							val.fo = realCryst;
							val.fc = atomReal;
							val.weight = atomImag;
							#ifdef COORDVAL_FULL
							vec3 frac = fftCrystal->fracFromElement(cIndex);
							val.pos = frac;
							#endif

							vals->push_back(val);
						}
					}
					else
					{
						vals->at(count).fc = atomReal;
						vals->at(count).weight = atomImag;
						count++;
					}
				}
				else if (mapScoreType == MapScoreTypeNone)
				{
					/* Add the density to the real value of the crystal
					 * voxel. */

					fftCrystal->addToReal(cIndex, atomReal);
				}
			}
		}
	}

	return 0;
}

void VagFFT::setScale(double cubeDim)
{
	mat3x3 real = make_mat3x3();
	mat3x3_scale(&real, cubeDim, cubeDim, cubeDim);
	_realBasis = real;
	adjustNs();

	_recipBasis = mat3x3_inverse(real);
	
	_toReal = _realBasis;
	mat3x3_scale(&_toReal, _nx, _ny, _nz);

	_toRecip = mat3x3_inverse(_toReal);
	_recipTrans = mat3x3_transpose(_toRecip);

	calculateAdjustmentVolumes();
	_setMatrices = true;
}

void VagFFT::addLessSimple(VagFFTPtr v2)
{
	for (int i = 0; i < nn(); i++)
	{
		long f1 = finalIndex(i);
		long f2 = v2->finalIndex(i);
		_data[f1][0] += v2->_data[f2][0];
	}

}

void VagFFT::addSimple(VagFFTPtr v2)
{
	if (v2->_total != _total)
	{
		addLessSimple(v2);
		return;
	}

	for (int i = 0; i < nn(); i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][0] += v2->_data[final_index][0];
	}
}

void VagFFT::copyRealToImaginary()
{
	for (int i = 0; i < _nn; i++)
	{
		long final_index = finalIndex(i);
		_data[final_index][1] = _data[final_index][0];
	}
}


double VagFFT::nanlessAverage()
{
	double sum = 0;
	double count = 0;

	for (int i = 0; i < _nn; i++)
	{
		long index = finalIndex(i);
		double val = _data[index][0];
		
		if (val != val)
		{
			continue;
		}

		sum += fabs(val);
		count++;
	}

	return sum / count;
}

double VagFFT::sumReal(int scratch)
{
	double sum = 0;

	for (int i = 0; i < _nn; i++)
	{
		long index = 0;
		if (scratch < 0 || scratch >= _nscratch)
		{
			index = finalIndex(i);
		}
		else
		{
			index = scratchIndex(i, scratch);
		}

		sum += fabs(_data[index][0]);
	}

	return sum;
}

double VagFFT::getIntensity(long i)
{
	long index = finalIndex(i);
	double val = _data[index][0] * _data[index][0] + 
	_data[index][1] * _data[index][1];
	return val;
}

double VagFFT::getAmplitude(long i)
{
	long index = finalIndex(i);
	double val = _data[index][0] * _data[index][0] + 
	_data[index][1] * _data[index][1];
	return sqrt(val);
}

vec3 VagFFT::fracFromElement(long int element)
{
	long x = element % _nx;
	element -= x;
	element /= _nx;

	long y = element % _ny;
	element -= y;
	element /= _ny;

	long z = element;

	double xfrac = (double)x / (double)_nx;
	double yfrac = (double)y / (double)_ny;
	double zfrac = (double)z / (double)_nz;

	return make_vec3(xfrac, yfrac, zfrac);
}

void VagFFT::writeToFile(std::string filename, double maxResolution,
                         VagFFTPtr data, VagFFTPtr diff, VagFFTPtr calc)
{
	if (_status != FFTReciprocalSpace)
	{
		std::cout << "Warning: writing MTZ file for grid in "
		"real space: " << filename << std::endl;
	}
	
	double nLimit[3];
	nLimit[0] = _nx;
	nLimit[1] = _ny;
	nLimit[2] = _nz;
	
	for (int i = 0; i < 3; i++)
	{
		nLimit[i] = nLimit[i] - ((int)nLimit[i] % 2); // make even
		nLimit[i] /= 2;
	}
	
	if (data)
	{
		// lower limit if data has smaller spacing
		nLimit[0] = (nLimit[0] > data->nx()) ? data->nx() : nLimit[0];
		nLimit[1] = (nLimit[1] > data->ny()) ? data->ny() : nLimit[1];
		nLimit[2] = (nLimit[2] > data->nz()) ? data->nz() : nLimit[2];
	}
	
	double dStar = 1 / maxResolution;

	if (maxResolution <= 0) dStar = FLT_MAX;

	/* For writing MTZ files */

	int columns = 12;

	float cell[6], wavelength;
	float *fdata = new float[columns];

	/* variables for symmetry */
	float rsm[192][4][4];
	char ltypex[2];

	/* variables for MTZ data structure */
	CMtz::MTZ *mtzout;
	CMtz::MTZXTAL *xtal;
	CMtz::MTZSET *set;
	CMtz::MTZCOL *colout[columns + 1];

	double unitCell[6];
	unit_cell_from_mat3x3(_toReal, unitCell);
	cell[0] = unitCell[0];
	cell[1] = unitCell[1];
	cell[2] = unitCell[2];
	cell[3] = unitCell[3];
	cell[4] = unitCell[4];
	cell[5] = unitCell[5];
	wavelength = 1.00; // fixme

	std::string outputFileOnly = filename;
	std::string outputFile = FileReader::addOutputDirectory(outputFileOnly);

	mtzout = CMtz::MtzMalloc(0, 0);
	ccp4_lwtitl(mtzout, "Written from Helen's XFEL tasks ", 0);
	mtzout->refs_in_memory = 0;
	mtzout->fileout = CMtz::MtzOpenForWrite(outputFile.c_str());

	if (!mtzout->fileout)
	{
		std::cout << "Could not open " << outputFile.c_str() << std::endl;
		std::cerr << "Error: " << strerror(errno) << std::endl;
		return;
	}

	// then add symm headers...
	for (int i = 0; i < _spg->nsymop; ++i)
	CCP4::rotandtrn_to_mat4(rsm[i], _spg->symop[i]);
	strncpy(ltypex, _spg->symbol_old, 1);
	ccp4_lwsymm(mtzout, _spg->nsymop, _spg->nsymop_prim, rsm, ltypex,
	            _spg->spg_ccp4_num, _spg->symbol_old, _spg->point_group);

	// then add xtals, datasets, cols
	xtal = MtzAddXtal(mtzout, "vagabond_crystal", "vagabond_project", cell);
	set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
	colout[0] = MtzAddColumn(mtzout, set, "H", "H");
	colout[1] = MtzAddColumn(mtzout, set, "K", "H");
	colout[2] = MtzAddColumn(mtzout, set, "L", "H");
	colout[3] = MtzAddColumn(mtzout, set, "FREE", "I");
	colout[4] = MtzAddColumn(mtzout, set, "FP", "F");
	colout[5] = MtzAddColumn(mtzout, set, "SIGFP", "Q");
	
	if (&*data != this)
	{
		colout[6] = MtzAddColumn(mtzout, set, "FC", "F");
		colout[7] = MtzAddColumn(mtzout, set, "FWT", "F");
		colout[8] = MtzAddColumn(mtzout, set, "PHIC", "P");
		colout[9] = MtzAddColumn(mtzout, set, "PHWT", "P");
		colout[10] = MtzAddColumn(mtzout, set, "DELFWT", "F");
		colout[11] = MtzAddColumn(mtzout, set, "PHDELWT", "P");
	}

	int num = 0;
	
	mat3x3 transpose = mat3x3_transpose(_toRecip);

	/* symmetry issues */
	for (int k = -nLimit[2]; k < nLimit[2]; k++)
	{
		for (int j = -nLimit[1]; j < nLimit[1]; j++)
		{
			for (int i = -nLimit[0]; i < nLimit[0]; i++)
			{
				bool asu = CSym::ccp4spg_is_in_asu(_spg, i, j, k);
				bool f000 = (i == 0 && j == 0 && k == 0);
				int abs = CSym::ccp4spg_is_sysabs(_spg, i, j, k);
				
				if (!asu || f000 || abs)
				{
					continue;
				}

				vec3 pos = make_vec3(i, j, k);
				mat3x3_mult_vec(transpose, &pos);

				if (vec3_length(pos) > dStar)
				{
					continue;
				}

				/* weighted amplitude and phase */
				double fwt = getAmplitude(i, j, k);
				double phwt = getPhase(i, j, k);

				/* get calculated amplitude from calc map, ideally
				 * available */
				double calcAmp = getAmplitude(i, j, k);
				double phic = phwt;
				
				if (calc)
				{
					calcAmp = calc->getAmplitude(i, j, k);
					phic = calc->getPhase(i, j, k);
				}

				double free = 1;

				/* Sort out observed values from observed diffraction,
				 * ideally available */
				double foInt = fwt * fwt;
				double foAmp = fwt;
				double sigma = 1;

				if (data)
				{
					int ele = data->element(i, j, k);
					foAmp = data->getReal(ele);
					sigma = data->getImag(ele);
					free = data->getScratchComponent(ele, 0, 0);
				}


				double diffPhwt = 0;
				double diffAmp = 0;
				
				if (diff)
				{
					diffPhwt = diff->getPhase(i, j, k);
					diffAmp = diff->getAmplitude(i, j, k);
				}

				/* MTZ file stuff */

				if (f000)
				{
					fwt = calcAmp;
					diffAmp = 0;
				}

				fdata[0] = i;
				fdata[1] = j;
				fdata[2] = k;
				fdata[3] = free - 0.1;
				fdata[4] = foAmp;
				fdata[5] = sigma;

				if (&*data != this)
				{
					fdata[6] = calcAmp;
					fdata[7] = fwt;
					fdata[8] = phic;
					fdata[9] = phwt;
					fdata[10] = diffAmp;
					fdata[11] = diffPhwt;
				}

				num++;
				ccp4_lwrefl(mtzout, fdata, colout, columns, num);
			}
		}
	}

	MtzPut(mtzout, " ");
	MtzFree(mtzout);

	delete [] fdata;

}

void VagFFT::applySymmetry(bool silent, double topRes, bool average)
{
	fftwf_complex *tempData;
	tempData = (fftwf_complex *)fftwf_malloc(_nn * sizeof(FFTW_DATA_TYPE));
	memset(tempData, 0, sizeof(FFTW_DATA_TYPE) * _nn);

	int count = 0;
	
	if (_spg->spg_num == 1)
	{
		return;
	}

	if (!silent)
	{
		std::cout << "applying " <<
		_spg->nsymop << " symmetry operators, space group " 
		<< _spg->symbol_xHM;
		std::cout << " (" << _spg->spg_num << ")"  << ": " << std::flush;
	}

	/* Loop through and convert data into amplitude and phase */
	for (int pre_n = 0; pre_n < _nn; pre_n++)
	{
		long n = finalIndex(pre_n);
		double xOrig = _data[n][0];
		double yOrig = _data[n][1];
		double myAmp = sqrt(xOrig * xOrig + yOrig * yOrig);
		double myPhase = atan2(yOrig, xOrig) * 180 / M_PI;
		while (myPhase >= 360) myPhase -= 360;
		while (myPhase < 0) myPhase += 360;

		_data[n][0] = myAmp;
		_data[n][1] = myPhase;
	}
	
	double dMax = FLT_MAX;
	if (topRes > 0)
	{
		dMax = 1 / topRes;
	}

	for (int k = -_nz / 2; k < _nz / 2; k++)
	{
		for (int j = -_ny / 2; j < _ny / 2; j++)
		{
			for (int i = -_nx / 2; i < _nx / 2; i++)
			{
				int abs = CSym::ccp4spg_is_sysabs(_spg, i, j, k);
				long pre_index = element(i, j, k);
				int sym = 0;
				
				if (!average)
				{
					int _h, _k, _l;
					sym = CSym::ccp4spg_put_in_asu(_spg, i, j, k, 
					                                   &_h, &_k, &_l);
					pre_index = element(_h, _k, _l);
				}

				long index = finalIndex(pre_index);

				if (abs)
				{
					/* we know this has to be zero, tempData already is */
					continue;	
				}
				
				if (topRes > 0)
				{
					vec3 ijk = make_vec3(i, j, k);
					mat3x3_mult_vec(_toRecip, &ijk);

					double sqlength = vec3_sqlength(ijk);

					if (sqlength > dMax * dMax)
					{
						continue;
					}
				}

				/* Not misnomers: dealt with in previous paragraph */
				double myAmp = _data[index][0];
				double myPhase = _data[index][1];

				/* if averaging, loop through all symops 
				 * this index applies to */

				for (int l = 0; l < _spg->nsymop && average; l++)
				{
					float *rot = &_spg->invsymop[l].rot[0][0];

					/* rotation */
					int _h, _k, _l;
					_h = (int) rint(i*rot[0] + j*rot[3] + k*rot[6]);
					_k = (int) rint(i*rot[1] + j*rot[4] + k*rot[7]);
					_l = (int) rint(i*rot[2] + j*rot[5] + k*rot[8]);

					long sym_index = element(_h, _k, _l);
					/* translation */
					float *trn = _spg->symop[l].trn;

					double shift = (float)_h * trn[0];
					shift += (float)_k * trn[1];
					shift += (float)_l * trn[2];
					shift = shift - floor(shift);

					double deg = myPhase + shift * 360.;
					double newPhase = deg2rad(deg);
					
					/* add to temporary data array */
					double x = myAmp * cos(newPhase);
					double y = myAmp * sin(newPhase);

					tempData[sym_index][0] += x;
					tempData[sym_index][1] += y;
				}
				
				/* otherwise, grab symops and apply to oneself */
				if (!average)
				{
					int symop = (sym - 1) / 2;

					float *trn = _spg->symop[symop].trn;

					/* translation */
					double shift = (float)i * trn[0];
					shift += (float)j * trn[1];
					shift += (float)k * trn[2];
					shift = shift - floor(shift);
					
					double deg = myPhase + shift * 360.;
					double newPhase = deg2rad(deg);
					
					/* add to temporary data array */
					double x = myAmp * cos(newPhase);
					double y = myAmp * sin(newPhase);

					long tmp_index = element(i, j, k);
					tempData[tmp_index][0] = myAmp;
					tempData[tmp_index][1] = 0;
				}
			}
		}
	}

	if (!silent)
	{
		std::cout << "... done." << std::endl;
	}

	/* Loop through and convert data into amplitude and phase */
	for (int n = 0; n < _nn; n++)
	{
		long m = finalIndex(n);
		_data[m][0] = tempData[n][0];
		_data[m][1] = tempData[n][1];
	}

	free(tempData);

}

void VagFFT::multiplyAll(float value)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		_data[index][0] *= value;
		_data[index][1] *= value;
	}
}

void VagFFT::copyFrom(VagFFTPtr other)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		int iother = other->finalIndex(i);
		_data[index][0] = other->_data[iother][0];
		_data[index][1] = other->_data[iother][1];
	}
}

void VagFFT::addToScratch(int scratch)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		int sindex = scratchIndex(i, scratch);
		_data[sindex][0] += _data[index][0];
		_data[sindex][1] += _data[index][1];
		_data[index][0] = 0;
		_data[index][1] = 0;
	}
}

void VagFFT::copyToScratch(int scratch)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		int sindex = scratchIndex(i, scratch);
		_data[sindex][0] = _data[index][0];
		_data[sindex][1] = _data[index][1];
		_data[index][0] = 0;
		_data[index][1] = 0;
	}
}

void VagFFT::addScratchBack(int scratch)
{
	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		int sindex = scratchIndex(i, scratch);
		_data[index][0] += _data[sindex][0];
		_data[index][1] += _data[sindex][1];

	}
}

double VagFFT::minValue()
{
	double min = FLT_MAX;

	for (long i = 0; i < nn(); i++)
	{
		int index = finalIndex(i);
		if (min > _data[index][0])
		{
			min = _data[index][0];
		}
	}
	
	return min;
}

double VagFFT::compareReciprocalToScratch(int scratch)
{
	double sum = 0;
	double stdev = 5.0;
	double weights = 0;

	for (int k = -_nz / 2; k < _nz / 2; k++)
	{
		for (int j = -_ny / 2; j < _ny / 2; j++)
		{
			for (int i = -_nx / 2; i < _nx / 2; i++)
			{
				int abs = CSym::ccp4spg_is_sysabs(_spg, i, j, k);

				if (abs)
				{
					continue;	
				}

				bool asu = CSym::ccp4spg_is_in_asu(_spg, i, j, k);
				
				if (!asu)
				{
					continue;
				}

				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(_toRecip, &ijk);

				double sqlength = vec3_sqlength(ijk);

				long pre_index = element(i, j, k);
				long scr_index = scratchIndex(pre_index, scratch);
				
				float y = _data[scr_index][0]; /* reference */
				
				if (y != y)
				{
					continue;
				}

				double x = getIntensity(pre_index);
				
				double weight = 1;

				double err = fabs(x - y) / y;
				err /= stdev;
				double e = exp(-(err * err));
				
				if (e != e)
				{
					continue;
				}
				
				sum += e;
				weights += weight;
			}
		}
	}
	
	/*
	std::cout << "Compare recip to scratch: " << 
	std::setprecision(6) << sum / weights << std::endl;
	*/

	return sum / weights;
}

void VagFFT::setAllReal(double val)
{
	for (long i = 0; i < nn(); i++)
	{
		long ii = finalIndex(i);
		_data[ii][0] = val;
		_data[ii][1] = 0;
	}
}

double VagFFT::sumAmp(int scratch)
{
	double sum = 0;

	for (long i = 0; i < nn(); i++)
	{
		int in = 0;
		if (scratch < 0)
		{
			in = finalIndex(i);
		}
		else
		{
			in = scratchIndex(i, scratch);
		}

		double intensity = _data[in][0] * _data[in][0] 
		+ _data[in][1] * _data[in][1];
		
		if (intensity != intensity)
		{
			continue;
		}
		
		sum += intensity;
	}
	
	return sum;
}

void VagFFT::findLimitingValues(double xMin, double xMax, double yMin,
                                double yMax, double zMin, double zMax,
                                vec3 *minVals, vec3 *maxVals)
{
	mat3x3 toCrystBasis = getRecipBasis();
	*minVals = make_vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	*maxVals = make_vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);

	for (double k = zMin; k <= zMax + 0.1; k += (zMax - zMin))
	{
		for (double j = yMin; j <= yMax + 0.1; j += (zMax - zMin))
		{
			for (double i = xMin; i <= xMax + 0.1; i += (zMax - zMin))
			{
				vec3 angCorner = make_vec3(i, j, k);
				vec3 toCrystal = mat3x3_mult_vec(toCrystBasis, angCorner);
				toCrystal.x = lrint(toCrystal.x);
				toCrystal.y = lrint(toCrystal.y);
				toCrystal.z = lrint(toCrystal.z);

				vec3_min_each(minVals, toCrystal);
				vec3_max_each(maxVals, toCrystal);
			}
		}
	}
}

/* is a given ijk within -n/2 < ijk <= n/2 */
bool VagFFT::withinBounds(int i, int j, int k)
{
	if (i >= _nx / 2 || i < -_nx / 2)
	{
		return false;
	}

	if (j >= _ny / 2 || j < -_ny / 2)
	{
		return false;
	}

	if (k >= _nz / 2 || k < -_nz / 2)
	{
		return false;
	}
	
	return true;
}


void VagFFT::setScratchComponent(int index, int scratch, int comp, double val)
{
	long si = scratchIndex(index, scratch);
	_data[si][comp] = val;
}

void VagFFT::shrink(double radius)
{
	vec3 mins = make_vec3(0, 0, 0);
	vec3 maxs = make_vec3(0, 0, 0);
	mat3x3 basis = getRealBasis();
	findLimitingValues(0, radius, 0, radius, 0, radius, &mins, &maxs);
					
	int count = 0;
	int total = 0;
	int solv = 0;

	for (int z = 0; z < nz(); z++)
	{
		for (int y = 0; y < ny(); y++)
		{
			for (int x = 0; x < nx(); x++)
			{
				long raw = element(x, y, z);
				long iraw = finalIndex(raw);
				/* We only want to modify protein to become
				 * more like solvent */
				if (_data[iraw][0] == 1)
				{
					solv++;
					_data[iraw][1] = 1;
					continue;
				}

				/* Default is to be protein */
				_data[raw][1] = 0;
				bool done = false;

				for (int k = 0; k < maxs.z && !done; k++)
				{
					for (int j = 0; j < maxs.y && !done; j++)
					{
						for (int i = 0; i < maxs.x && !done; i++)
						{
							vec3 ijk = make_vec3(i, j, k);
							mat3x3_mult_vec(basis, &ijk);

							/* Doesn't matter if it goes over radial
							 * boundary */
							if (vec3_sqlength(ijk) > radius * radius)
							{
								continue;
							}

							long index = element(i + x, j + y, k + z);
							long ii = finalIndex(index);
							
							if (_data[ii][0] == 1)
							{
								count++;
								done = true;
								_data[iraw][1] = 1;
								break;
							}
						}
					}
				}
				
				total++;
			}
		}
	}
	
	for (int i = 0; i < nn(); i++)
	{
		_data[i][0] = _data[i][1];
		_data[i][1] = 0;
	}
}


void VagFFT::bittyShrink(double radius, int num)
{
	vec3 mins = make_vec3(0, 0, 0);
	vec3 maxs = make_vec3(0, 0, 0);
	mat3x3 basis = getRealBasis();
	findLimitingValues(0, radius, 0, radius, 0, radius,
					   &mins, &maxs);

	for (int z = 0; z < nz(); z++)
	{
		for (int y = 0; y < ny(); y++)
		{
			for (int x = 0; x < nx(); x++)
			{
				long raw = element(x, y, z);
				long iraw = finalIndex(raw);
				bool change = false;

				int work = (int)_data[iraw][0];
				
				for (int b = 0; b < num; b++)
				{
					unsigned char byte = (work >> b) & 1;
					unsigned char newbyte = 1; /* protein */

					/* We only want to modify protein to become
					 * more like solvent, so ignore solvent */
					if (byte == 0)
					{
						/* set new bit to solvent and continue */
						work &= ~(1 << b + num); 
						continue;
					}

					bool done = false;

					/* Default is to be protein */
					newbyte = 1;

					/* now byte is protein, but will it remain protein? */
					for (int k = 0; k < maxs.z && !done; k++)
					{
						for (int j = 0; j < maxs.y && !done; j++)
						{
							for (int i = 0; i < maxs.x && !done; i++)
							{
								vec3 ijk = make_vec3(i, j, k);

								mat3x3_mult_vec(basis, &ijk);

								/* Doesn't matter if it goes over radial
								 * boundary */
								if (vec3_sqlength(ijk) > radius * radius)
								{
									continue;
								}

								long index = element(i + x, j + y, k + z);
								long ii = finalIndex(index);
								int neigh = _data[ii][0];

								/* Switch to solvent if we found solvent */
								unsigned char other = (neigh >> b) & 1;

								if (other == 0)
								{
									done = true;
									change = true;
									newbyte = 0;
									break;
								}
							}
						}
					}
					
					work |= (newbyte << b + num);
				}
				
				int newmask = work >> num;
				_data[iraw][1] = (float)newmask;
			}
		}
	}
	
	for (int i = 0; i < nn(); i++)
	{
		_data[i][0] = _data[i][1];
		_data[i][1] = 0;
	}
}

void printb(int val)
{
	for (int i = 16; i >= 0; i--)
	{
		std::cout << ((val >> i) & 0b00000001 ? 1 : 0);
	}

	std::cout << std::endl;
}

void VagFFT::addToValueAroundPoint(vec3 pos, double radius, double value,
                                   int bitIndex)
{
	/* Determine square bounding box for radius in Ang. */
	vec3 minRadius, maxRadius;
	
	mat3x3_mult_vec(_toRecip, &pos);
	collapseFrac(&pos.x, &pos.y, &pos.z);
	pos.x *= nx();
	pos.y *= ny();
	pos.z *= nz();
	
	findLimitingValues(-radius, radius, -radius, radius, -radius, radius,
	                   &minRadius, &maxRadius);
	mat3x3 basis = getRealBasis();
	
	bool useBit = (bitIndex >= 0);
	int bitty = 0;
	if (useBit)
	{
		bitty = 1 << bitIndex;
	}
	
	for (double k = minRadius.z; k < maxRadius.z; k++)
	{
		for (double j = minRadius.y; j < maxRadius.y; j++)
		{
			for (double i = minRadius.x; i < maxRadius.x; i++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(basis, &ijk);
				
				/* Doesn't matter if it goes over radial
				 * boundary */
				if (vec3_sqlength(ijk) > radius * radius)
				{
					continue;
				}

				double x = lrint(i + pos.x);
				double y = lrint(j + pos.y);
				double z = lrint(k + pos.z);
				
				long index = element(x, y, z);
				long ii = finalIndex(index);

				if (!useBit)
				{
					_data[ii][0] += value;

					if (_data[ii][0] < 0)
					{
						_data[ii][0] = 0;
					}
				}
				else
				{
					int work = (int)_data[ii][0];
					work |= bitty;
					
					_data[ii][0] = (float)work;
				}
			}
		}
	}
}

vec3 VagFFT::getSymRelatedPosition(vec3 pos, int i)
{
	mat3x3_mult_vec(_toRecip, &pos);

	float *rot = &_spg->symop[i].rot[0][0];
	float *trn = _spg->symop[i].trn;

	vec3 mod = empty_vec3();
	mod.x = pos.x * rot[0] + pos.y * rot[1] + pos.z * rot[2];
	mod.y = pos.x * rot[3] + pos.y * rot[4] + pos.z * rot[5];
	mod.z = pos.x * rot[6] + pos.y * rot[7] + pos.z * rot[8];
	mod.x += trn[0];
	mod.y += trn[1];
	mod.z += trn[2];
	
	mat3x3_mult_vec(_toReal, &mod);
	
	return mod;
}

void VagFFT::convertMaskToSolvent(int expTotal)
{
	for (size_t i = 0; i < nn(); i++)
	{
		int ayes = 0;
		int noes = 0;

		long ii = finalIndex(i);
		int val = _data[ii][0];
		
		for (int j = 0; j < expTotal; j++)
		{
			unsigned char byte = (val >> j) & 1;
			
			if (byte)
			{
				ayes++;
			}
			else
			{
				noes++;
			}
		}

		/* no is solvent */
		float frac = (float)noes / (float)(noes + ayes);
		_data[ii][0] = frac;
		_data[ii][1] = 0; /* should be the case */
	}
}

void VagFFT::reindex(mat3x3 reindex)
{
	fftwf_complex *tempData;
	tempData = (fftwf_complex *)fftwf_malloc(_nn * sizeof(FFTW_DATA_TYPE));
	memset(tempData, 0, sizeof(FFTW_DATA_TYPE) * _nn);

	for (int k = -_nz / 2; k < _nz / 2; k++)
	{
		for (int j = -_ny / 2; j < _ny / 2; j++)
		{
			for (int i = -_nx / 2; i < _nx / 2; i++)
			{
				vec3 ijk = make_vec3(i, j, k);
				mat3x3_mult_vec(reindex, &ijk);
				
				if (!withinBounds(ijk.x, ijk.y, ijk.z))
				{
					continue;
				}
				
				long end = element(ijk.x, ijk.y, ijk.z);

				long ijk_index = element(i, j, k);
				long start = finalIndex(ijk_index);
				
				tempData[end][0] = _data[start][0];
				tempData[end][1] = _data[start][1];
			}
		}
	}

	/* Loop through and convert data into amplitude and phase */
	for (int n = 0; n < _nn; n++)
	{
		long m = finalIndex(n);
		_data[m][0] = tempData[n][0];
		_data[m][1] = tempData[n][1];
	}

	free(tempData);
}

double VagFFT::resolution(int i, int j, int k)
{
	vec3 ijk = make_vec3(i, j, k);
	mat3x3_mult_vec(_recipTrans, &ijk);
	return 1 / vec3_length(ijk);
}

double VagFFT::voxelVolume()
{
	return mat3x3_volume(_realBasis);
}
