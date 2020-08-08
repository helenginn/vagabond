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

class AveCSV : public Average
{
public:
	AveCSV(Group *group, std::string csv);

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
	void load();
	virtual void calculate();
	virtual double findCorrelation(MtzFFTPtr one, MtzFFTPtr two);
	void addValue(std::string id0, std::string id1, double val);
private:
	std::string _csv;
	static bool _usingCSV;

	ClusterList *_list;
	std::map<std::string, int> _ids;
	static std::map<std::string, MapStringDouble> _relationships;
};

#endif
