// 
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

#include "Balance.h"
#include "Bond.h"
#include "ParamBand.h"

Balance::Balance(BondPtr bond)
{
	if (bond->downstreamBondGroupCount() <= 1)
	{
		return;
	}
	
	_bond = bond;
	
	for (int i = 0; i < bond->downstreamBondGroupCount(); i++)
	{
		ParamBandPtr band = ParamBandPtr(new ParamBand());
		band->setPrivateGetter(Bond::getOccupancy);
		band->setPrivateSetter(Bond::setOccupancy);
		band->setMinMax(0.05, 1.);

		/* Don't want to divide occupancy out between multiple bond groups */
		band->setEqualDivision(false);

		for (int j = 0; j < bond->downstreamBondCount(i); j++)
		{
			BondPtr tmp = bond->downstreamBond(i, j);
			
			if (tmp->isSplit())
			{
				band->addObject(&*bond->downstreamBond(i, j), 1);
			}
		}
		
		band->prepare();
		
		if (band->objectCount())
		{
			_bands.push_back(band);
		}
	}
}

double Balance::addParamsToStrategy(RefinementStrategyPtr strategy)
{
	for (int i = 0; i < _bands.size(); i++)
	{
		strategy->addParameter(&*_bands[i], ParamBand::getGlobalParam,
		                       ParamBand::setGlobalParam, 0.2, 0.01);
	}
}

void Balance::adjustment()
{
	double sum = 0;

	for (int i = 0; i < _bands.size(); i++)
	{
		ParamBandPtr band = _bands[i];
		double val = ParamBand::getGlobalParam(&*band) + 
		band->baseValueForObject(0);
		sum += val;
	}
	
	for (int i = 0; i < _bands.size(); i++)
	{
		ParamBandPtr band = _bands[i];
		double val = band->getParamValue();
		band->setWeightForAll(1 / sum);
		ParamBand::setGlobalParam(&*band, val);
	}
	
	_bond->propagateChange(10);
}

