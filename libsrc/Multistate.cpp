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

#include "Multistate.h"
#include <hcsrc/FileReader.h>
#include "libsrc/Crystal.h"
#include "libsrc/PDBReader.h"

Multistate::Multistate(std::string filename)
{
	_filename = filename;
}

void Multistate::process()
{
	if (!file_exists(_filename))
	{
		std::cout << "File does not exist: " << _filename << std::endl;
		return;
	}

	std::string contents = get_file_contents(_filename);
	std::stringstream ss(contents);

	std::string line;
	bool toModel = false;
	std::string model;
	while (std::getline(ss, line))
	{
		if (line.substr(0, 6) == "CRYST1")
		{
			_cryst1 = line;
		}
		else if (line.substr(0, 6) == "MODEL ")
		{
			toModel = true;
			model += _cryst1;
			model += "\n";
		}
		
		else if (line.substr(0, 6) == "ENDMDL")
		{
			if (model.length())
			{
				_individuals.push_back(model);
			}
			
			model = "";
			toModel = false;
		}
		else if (toModel)
		{
			model += line;
			model += "\n";
		}
	}
	
	for (size_t i = 0; i < _individuals.size(); i++)
	{
		PDBReader reader;
		reader.setContents(_individuals[i]);
		if (_atom.length())
		{
			reader.ignoreAtomsExcept(_atom);
		}

		CrystalPtr crystal = reader.getCrystal();
		if (crystal)
		{
			_crystals.push_back(crystal);
		}
	}
	
	std::cout << "Multistate PDB with " <<  _crystals.size() 
	<< " models." << std::endl;
}
