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

#ifndef __cluster4x__AveVectors__
#define __cluster4x__AveVectors__

#include "Average.h"

#include <string>
#include <map>

class ClusterList;

class AveVectors : public Average
{
public:
	AveVectors(Group *group, std::string csv = "");
	
	void setFilename(std::string file)
	{
		_csv = file;
	}
	
	void setList(ClusterList *list)
	{
		_list = list;
	}
	
	static bool hasVectorData()
	{
		return _vectors.size() > 0;
	}
	
	size_t titleCount()
	{
		return _titles.size();
	}
	
	void addTitle(std::string title)
	{
		_titles.push_back(title);
	}
	
	std::string title(int i)
	{
		return _titles[i];
	}
	
	size_t vectorCount()
	{
		return _vectors.size();
	}
	
	void setEnabled(int i, bool e)
	{
		_enabled[i] = e;
	}
	
	bool enabled(int i)
	{
		return _enabled[i];
	}
	
	std::vector<double> &vector(int i)
	{
		return _vectors[_names[i]];
	}

	void load();
	void preparePaths();
	void setVector(std::string name, std::vector<double> vec);
	virtual double findCorrelation(MtzFFTPtr one, MtzFFTPtr two);
	virtual void calculate();
private:
	void loadTitles(std::string line);
	std::string _csv;

	static std::vector<std::string> _names;
	static std::vector<std::string> _titles;
	static std::map<int, bool> _enabled;
	static std::map<std::string, int> _ids;
	static std::map<std::string, std::vector<double> > _vectors;
	std::vector<double> _averageVec;
	std::vector<double> _sigVec;
	ClusterList *_list;
};

#endif
