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

#ifndef __cluster4x__AveCSV__
#define __cluster4x__AveCSV__

#include <string>
#include <map>
#include "Average.h"

class ClusterList;

typedef std::map<std::string, double> MapStringDouble;
typedef std::map<std::string, MapStringDouble> Relationships;

class AveCSV : public Average
{
public:
	AveCSV(Group *group, std::string csv = "");
	
	void setFilename(std::string file)
	{
		_csv = file;
	}

	static bool usingCSV()	
	{
		return _usingCSV;
	}
	
	static void setUsingCSV(bool csv)
	{
		_usingCSV = csv;
	}

	void setList(ClusterList *list)
	{
		_list = list;
	}
	
	static void clear();
	
	static size_t csvCount()
	{
		return _relationships.size();
	}
	
	static std::string csvName(int i)
	{
		return _filenames[i];
	}
	
	static int currentChoice()
	{
		return _chosen;
	}
	
	static void setChosen(int chosen)
	{
		_chosen = chosen;
	}
	
	void setComparisonType(int type)
	{
		_compType = type;
	}
	
	static void setChosen(std::string);
	
	void startNewCSV(std::string name);

	void load();
	virtual void calculate();
	virtual double findCorrelation(MtzFFTPtr one, MtzFFTPtr two);
	void addValue(std::string id0, std::string id1, double val);
	void preparePaths();
private:
	std::string _csv;
	static bool _usingCSV;

	ClusterList *_list;
	std::map<std::string, int> _ids;
	std::map<std::string, std::vector<double> > _vectors;
	static std::vector<std::string> _filenames;
	static std::vector<Relationships> _relationships;
	static int _chosen;
	int _compType;
};

#endif
