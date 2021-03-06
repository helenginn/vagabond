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
#include "Options.h"
#include "Anchor.h"
#include "shared_ptrs.h"

void BucketPerStrand::addSolvent()
{
	CrystalPtr crystal = getCrystal();
	
	/* We run without a conformer to get our atom pointers from
	 * average values */
	addSolventForConformer(-1);

	VagFFTPtr f = crystal->getFFT();
	_solvent = VagFFTPtr(new VagFFT(*f, 1));
	_solvent->makePlans();
	_solvent->setStatus(FFTRealSpace);

	int confs = Anchor::sampleSize();
	double count = 0;

	std::cout << "Adding solvent for conformer " << std::flush;

	Timer t = Timer("adding solvent for conformers", true);
	for (int i = 0; i < crystal->getSampleNum(); i += SOLVENT_BITS)
	{
		std::cout << i << " ... " << std::flush;
		int num = SOLVENT_BITS;

		if (i + num > confs)
		{
			num = confs - i;
		}

		addSolventForConformer(i, num);
		/* does not appear to help matters here */
		
		double shrink = Options::getShrink();

		if (shrink > 0)
		{
			_solvent->bittyShrink(0.4, num);
		}
		
		_solvent->convertMaskToSolvent(num);
		_solvent->addToScratch(0);
		_solvent->multiplyAll(0);
		count++;
	}
	
	std::cout << std::endl;
	t.report();

	/* But overwrite the solvent afterwards with the new one we 
	 * have calculated */
	_solvent->addScratchBack(0);
	_solvent->multiplyAll(1 / (double)count);
	removeSlivers(1.5);
	
	reportSolventContent();
	adjustForVoxelVolume();
	
	_solvent->drawSlice(0, "solvent_slice_unfit");
	setPartialStructure(_solvent);
}
