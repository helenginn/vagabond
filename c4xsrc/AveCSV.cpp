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

#include "AveCSV.h"
#include "DatasetPath.h"
#include "ClusterList.h"
#include "MtzFile.h"
#include "MtzFFT.h"
#include <iostream>
#include <libsrc/FileReader.h>

std::map<std::string, MapStringDouble> AveCSV::_relationships;
bool AveCSV::_usingCSV = false;


AveCSV::AveCSV(Group *group, std::string csv) : Average(group)
{
	_symmetric = false;
	_list = NULL;
	_csv = csv;
	_usingCSV = true;
}

void AveCSV::load()
{
	std::string contents = get_file_contents(_csv);
	
	std::vector<std::string> lines = split(contents, '\n');
	
	for (size_t i = 0; i < lines.size(); i++)
	{
		std::cout << "Processing: " << lines[i] << std::endl;
		
		std::vector<std::string> components = split(lines[i], ',');
		
		if (components.size() != 3)
		{
			std::cout << "Wrong number of terms, skipping. " << std::endl;
			continue;
		}

		char *pos = NULL;
		double val = strtod(&components[2][0], &pos);

		if (pos == &(components[2][0]))
		{
			if (i == 0)
			{
				std::cout << "First line appears to be header." << std::endl;
			}
			else
			{
				std::cout << "Final value of line " << i << " is not "
				"a number?" << std::endl;
			}

			continue;
		}
		
		val = 1 - val;
		
		for (size_t j = 0; j < components.size(); j++)
		{
			std::string id0 = components[0];
			_ids[id0]++;
			std::string id1 = components[1];
			_ids[id1]++;
			
			_relationships[id0][id1] = val;
		}
	}
	
	std::vector<DatasetPath> paths;
	for (std::map<std::string, int>::iterator it = _ids.begin(); 
	     it != _ids.end(); it++)
	{
		std::string id = it->first;
		DatasetPath path;
		path.mtz_path = "";
		path.pandda_mtz = "";
		path.pdb_path = "";
		path.metadata = id;
		path.refinement_id = id;
		paths.push_back(path);
		std::cout << "Made entry for " << path.metadata << std::endl;
	}

	_list->load(paths);
}

void AveCSV::calculate()
{

}

double AveCSV::findCorrelation(MtzFFTPtr one, MtzFFTPtr two)
{
	std::string met1 = one->getMtzFile()->metadata();
	std::string met2 = two->getMtzFile()->metadata();
	
	if (_relationships.count(met1))
	{
		return _relationships[met1][met2];
	}

	return 0;
}
