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

#include "ErroneousZone.h"
#include "Crystal.h"
#include "Polymer.h"
#include "Monomer.h"

ErroneousZone::ErroneousZone(QWidget *parent, CrystalPtr crystal) 
: QMainWindow(parent)
{
	QString title = "Erroneous zones";
	setWindowTitle(title);
	resize(600, 400);
	
	_crystal = crystal;
	_polymer = ToPolymerPtr(_crystal->molecule(0));

	std::cout << "Autosetting polymer to " << 
	_polymer->getChainID() << std::endl;
	
	_monomers.push_back(_polymer->getMonomer(189));
	_monomers.push_back(_polymer->getMonomer(190));
	_monomers.push_back(_polymer->getMonomer(191));
	_monomers.push_back(_polymer->getMonomer(192));
	_monomers.push_back(_polymer->getMonomer(225));
	_monomers.push_back(_polymer->getMonomer(226));
	_monomers.push_back(_polymer->getMonomer(227));
	_monomers.push_back(_polymer->getMonomer(228));
	
	_everything = AtomGroupPtr(new AtomGroup());
	
	for (int i = 0; i < _monomers.size(); i++)
	{
		_everything->addAtomsFrom(_monomers[i]);
	}
	
	std::cout << "Erroneous zone has " << _everything->atomCount()
	<< " atoms." << std::endl;
	
	_start2FoFc = _everything->scoreWithMap(ScoreTypeAddDensity, _crystal);
	std::cout << "Weighted density: " << _start2FoFc << std::endl;
	
	for (int i = 0; i < _monomers.size(); i++)
	{
		_monomers[i]->setWeighting(0);
	}
	omitResidues(0.7);

	Options::getActiveCrystal()->wrapUpRefinement();

	return;
	for (double weight = 1; weight >= 0; weight -= 0.1)
	{
		omitResidues(weight);
	}
}

void ErroneousZone::omitResidues(double weight)
{
	for (int i = 0; i < _monomers.size(); i++)
	{
		_monomers[i]->setWeighting(weight);
	}
	
	_crystal->silentConcludeRefinement();
	
	double neg_diff = _everything->scoreWithMap(ScoreTypeAddDensity, _crystal,
	                                            "", MapScoreFlagDifference
	                                            | MapScoreFlagNegOnly);
	 
	double pos_diff = _everything->scoreWithMap(ScoreTypeAddDensity, _crystal,
	                                            "", MapScoreFlagDifference
	                                            | MapScoreFlagPosOnly);

	std::cout << "Weight: " << weight << " " << neg_diff << " " 
	<< neg_diff / _start2FoFc << " "
	<< pos_diff << " " << pos_diff / _start2FoFc << std::endl;

}
