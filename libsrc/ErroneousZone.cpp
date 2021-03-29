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

#include "ErroneousZone.h"
#include "Polymer.h"
#include "Monomer.h"
#include "Backbone.h"
#include "Sidechain.h"

ErroneousZone::ErroneousZone(PolymerPtr pol, int start, int end)
{
	_polymer = pol;
	_start = start;
	_end = end;
}

void ErroneousZone::calcDifferenceDensities()
{
	CrystalPtr target = Options::getActiveCrystal();

	for (int i = _start; i < _end; i++)
	{
		MonomerPtr mon = _polymer->getMonomer(i);
		
		if (!mon)
		{
			continue;
		}

		BackbonePtr backbone = mon->getBackbone();
		SidechainPtr sidechain = mon->getSidechain();
		double cc = -sidechain->scoreWithMap(ScoreTypeCorrel, target);
	}

}
