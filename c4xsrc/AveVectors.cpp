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
				continue;
			}
			else
			{
				std::cout << "Second value of line " << i << " is not "
				"a number?" << std::endl;
			}
		}
		
		std::vector<double> vec;
		for (size_t j = 1; j < components.size(); j++)
		{
			double val = strtod(&components[j][0], &pos);
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
	if (MyDictator::valueForKey("transpose") == "true")
	{
		transpose();
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

void AveVectors::calculate()
{
	if (!group())
	{
		setGroup(Group::topGroup());
	}

	std::vector<double> counts;
	_averageVec.clear();
	_sigVec.clear();
	
	std::vector<std::string> currentNames;

	if (group())
	{
		for (size_t i = 0; i < group()->mtzCount(); i++)
		{
			MtzFFTPtr one = group()->getMtz(i);
			std::string name = one->getMtzFile()->metadata();

			currentNames.push_back(name);
		}
	}
	else
	{
		currentNames = _names;
	}

	for (size_t i = 0; i < currentNames.size(); i++)
	{
		std::string name = currentNames[i];

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

		if (MyDictator::valueForKey("stdev") != "true")
		{
			_sigVec[j] = 1;
		}

		if (MyDictator::valueForKey("average") == "false")
		{
			_averageVec[j] = 0;
		}
	}
}

std::vector<double> AveVectors::adjustedVector(int j)
{
	std::vector<double> v = vector(j);

	if (_averageVec.size() < v.size() || _sigVec.size() < v.size())
	{
		calculate();
	}
	
	for (size_t i = 0; i < v.size(); i++)
	{
		v[i] -= _averageVec[i];
		v[i] /= _sigVec[i];
	}

	return v;
}

void AveVectors::incorporateAverages()
{
	if (_averageVec.size() == 0)
	{
		calculate();
	}

	for (size_t i = 0; i < _names.size(); i++)
	{
		std::vector<double> vec = _vectors[_names[i]];
		
		for (size_t j = 0; j < _titles.size(); j++)
		{
			double mean = _averageVec[j];
			double stdev = _sigVec[j];
			if (MyDictator::valueForKey("average") == "false")
			{
				mean = 0;
			}

			if (MyDictator::valueForKey("stdev") != "true")
			{
				stdev = 1;
			}

			double x = (vec[j] - mean) / stdev;
			std::cout << vec[j] << " becomes " << x << std::endl;
			_vectors[_names[i]][j] = x;
		}
	}
}

double AveVectors::findCorrelation(MtzFFTPtr one, MtzFFTPtr two)
{
	if (_averageVec.size() == 0)
	{
		calculate();
	}

	std::string f1 = one->getMtzFile()->metadata();
	std::string f2 = two->getMtzFile()->metadata();
	
	std::vector<double> vec1 = _vectors[f1];
	std::vector<double> vec2 = _vectors[f2];

	CorrelData cd = empty_CD();
	for (size_t i = 0; i < vec1.size(); i++)
	{
		if (_enabled.count(i) && !_enabled[i])
		{
			continue;
		}

		double mean = _averageVec[i];
		double stdev = _sigVec[i];

		if (MyDictator::valueForKey("average") == "false")
		{
			mean = 0;
		}

		if (MyDictator::valueForKey("stdev") != "true")
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

	cd.sum_x = 0;
	cd.sum_y = 0;
	
	if (MyDictator::valueForKey("stdev") != "true")
	{
//		double top = cd.sum_w * cd.sum_xy - cd.sum_x * cd.sum_y;
//		return top;
	}

	return evaluate_CD(cd);
}

void AveVectors::clearVectors()
{
	_names.clear();
	_titles.clear();
	_enabled.clear();
	_ids.clear();
	_vectors.clear();
}

void AveVectors::transpose()
{
	_averageVec.clear();
	_sigVec.clear();
	incorporateAverages();
	std::vector<std::string> newNames, newTitles;
	std::map<std::string, int> newIds;
	std::map<int, bool> newEnabled;
	newNames = _titles;
	newTitles = _names;
	std::map<std::string, std::vector<double> > newVectors;

	for (size_t i = 0; i < _titles.size(); i++)
	{
		newIds[_titles[i]]++;
		
		std::vector<double> newVec;
		
		for (size_t j = 0; j < _names.size(); j++)
		{
			std::vector<double> oldVec = _vectors[_names[j]];
			newVec.push_back(oldVec[i]);
		}

		newVectors[_titles[i]] = newVec;
		newEnabled[i] = true;
	}
	
	_names = newNames;
	_titles = newTitles;
	_ids = newIds;
	_enabled = newEnabled;
	_vectors = newVectors;

	_averageVec.clear();
	_sigVec.clear();
}

void AveVectors::exportValues(std::string filename)	
{
	if (_averageVec.size() == 0)
	{
		calculate();
	}

	std::ofstream file;
	file.open(filename);
	
	file << ",";

	for (size_t j = 0; j < _titles.size(); j++)
	{
		file << _titles[j] << ",";
	}
	
	file << std::endl;

	for (size_t i = 0; i < _names.size(); i++)
	{
		std::vector<double> vec = _vectors[_names[i]];
		file << _names[i] << ",";
		
		for (size_t j = 0; j < _titles.size(); j++)
		{
			double mean = _averageVec[j];
			double stdev = _sigVec[j];
			if (MyDictator::valueForKey("average") == "false")
			{
				mean = 0;
			}

			if (MyDictator::valueForKey("stdev") != "true")
			{
				stdev = 1;
			}

			double x = (vec[j] - mean) / stdev;
			file << x << ",";
		}

		file << std::endl;
	}

	file.close();
}

void AveVectors::setAveDev(bool ave, bool dev)
{
	std::string aveStr = ave ? "true" : "false";
	std::string devStr = dev ? "true" : "false";
	MyDictator::setValueForKey("average", aveStr);
	MyDictator::setValueForKey("stdev", devStr);

}
