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

#include "BucketPerStrand.h"
#include "Crystal.h"
#include "shared_ptrs.h"
#include "fftw3d.h"

void BucketPerStrand::addSolvent()
{
	std::cout << "Attempting to add solvent." << std::endl;
	CrystalPtr crystal = getCrystal();

	FFTPtr total = FFTPtr(new FFT(*crystal->getFFT()));
	total->setAll(0);

	int confs = crystal->getSampleNum();
	int count = 0;

	for (int i = 0; i < crystal->getSampleNum(); i += SOLVENT_BITS)
	{
		int num = SOLVENT_BITS;

		if (i + num > confs)
		{
			num = confs - i;
		}

		std::cout << "Adding solvent for conformer " << i << "." << std::endl;
		addSolventForConformer(i, num);
		_solvent->bittyShrink(0.4, num);
		_solvent->convertMaskToSolvent(num);
		reportSolventContent();
		FFT::addSimple(total, _solvent);
		count++;
	}
	
	/* We run without a conformer to get our atom pointers from
	 * average values */
	addSolventForConformer(-1);

	/* But overwrite the solvent afterwards with the new one we 
	 * have calculated */
	_solvent = total;
	_solvent->multiplyAll(1 / (double)count);
	removeSlivers(2.0);
	
	reportSolventContent();
	
	setPartialStructure(_solvent);
}
