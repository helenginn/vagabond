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

#include <hcsrc/FileReader.h>
#include <hcsrc/maths.h>
#include "MyDictator.h"
#include "Group.h"
#include "MtzFFT.h"
#include "MtzFile.h"
#include "AveVectors.h"
#include "DatasetPath.h"
#include "ClusterList.h"

std::map<int, bool> AveVectors::_enabled;
std::vector<std::string> AveVectors::_names;
std::vector<std::string> AveVectors::_titles;
std::map<std::string, int> AveVectors::_ids;
std::map<std::string, std::vector<double> > AveVectors::_vectors;

AveVectors::AveVectors(Group *group, std::string csv) : Average(group)
{
	_csv = csv;
	_list = NULL;
	_vectorList = true;
}

void AveVectors::loadTitles(std::string line)
{
	std::vector<std::string> components = split(line, ',');
	
	for (size_t i = 1; i < components.size(); i++)
	{
		_titles.push_back(components[i]);
	}
}

void AveVectors::load()
{
	std::string contents = get_file_contents(_csv);
	std::vector<std::string> lines = split(contents, '\n');
	int max = 0;

	for (size_t i = 0; i < lines.size(); i++)
	{
		std::cout << "Processing: " << lines[i] << std::endl;

		std::vector<std::string> components = split(lines[i], ',');

		if (components.size() <= 1)
		{
			std::cout << "Wrong number of terms, skipping. " << std::endl;
			continue;
		}
		
		std::string name = components[0];

		char *pos = NULL;
		strtod(&components[1][0], &pos);

		if (pos == &(components[1][0]))
		{
			if (i == 0)
			{
				std::cout << "First line appears to be header." 
				<< std::endl;
				loadTitles(lines[0]);
			}
			else
			{
				std::cout << "Second value of line " << i << " is not "
				"a number?" << std::endl;
			}

			continue;
		}
		
		std::vector<double> vec;
		for (size_t j = 1; j < components.size(); j++)
		{
			double val = strtod(&components[j][0], &pos);
			std::cout << val << std::endl;
			vec.push_back(val);
			if (max < j)
			{
				max = j;
			}
			
			if (max > _titles.size())
			{
				_titles.resize(max);
			}
		}
		
		if (_ids[name] >= 1)
		{
			std::cout << "WARNING: duplicate sample " << name << std::endl;
			std::cout << "Overwriting!!" << std::endl;
			std::cout << std::endl;
		}
		
		setVector(name, vec);
	}
	
	for (size_t i = 0; i < _titles.size(); i++)
	{
		_enabled[i] = true;
	}
}

void AveVectors::setVector(std::string name, std::vector<double> vec)
{
	_ids[name]++;
	_vectors[name] = vec;
	_names.push_back(name);
}

void AveVectors::preparePaths()
{
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

void AveVectors::calculate()
{
	if (!group())
	{
		setGroup(Group::topGroup());
	}

	std::vector<double> counts;
	_averageVec.clear();
	_sigVec.clear();
	
	for (size_t i = 0; i < group()->mtzCount(); i++)
	{
		MtzFFTPtr one = group()->getMtz(i);
		std::string name = one->getMtzFile()->metadata();

		if (_vectors.count(name) == 0)
		{
			continue;
		}

		std::vector<double> vec = _vectors[name];
		
		if (vec.size() > _averageVec.size())
		{
			int start = _averageVec.size();
			_averageVec.resize(vec.size());
			_sigVec.resize(vec.size());
			counts.resize(vec.size());

			for (size_t j = start; j < _averageVec.size(); j++)
			{
				_sigVec[j] = 0;
				_averageVec[j] = 0;
				counts[j] = 0;
			}
		}

		for (size_t j = 0; j < vec.size(); j++)
		{
			double val = vec[j];
			if (val != val)
			{
				continue;
			}
			_averageVec[j] += val;
			_sigVec[j] += val * val;
			counts[j]++;
		}
	}
	
	for (size_t j = 0; j < _averageVec.size(); j++)
	{
		double sq = _sigVec[j];
		double sum = _averageVec[j];

		_sigVec[j] = sqrt(sq - sum * sum / counts[j]);
		_averageVec[j] /= counts[j];
	}
}

double AveVectors::findCorrelation(MtzFFTPtr one, MtzFFTPtr two)
{
	std::string f1 = one->getMtzFile()->metadata();
	std::string f2 = two->getMtzFile()->metadata();
	
	std::vector<double> vec1 = _vectors[f1];
	std::vector<double> vec2 = _vectors[f2];

	CorrelData cd = empty_CD();
	for (size_t i = 0; i < vec1.size(); i++)
	{
		if (!_enabled[i])
		{
			continue;
		}

		double mean = _averageVec[i];
		double stdev = _sigVec[i];

		if (MyDictator::valueForKey("average") == "false")
		{
			mean = 0;
		}

		if (MyDictator::valueForKey("stdev") == "false")
		{
			stdev = 1;
		}

		double x = (vec1[i] - mean) / stdev;
		double y = (vec2[i] - mean) / stdev;
		if (x != x || y != y)
		{
			continue;
		}

		add_to_CD(&cd, x, y);
	}

	return evaluate_CD(cd);
}
