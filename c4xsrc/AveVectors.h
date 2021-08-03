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

	void load();
	void preparePaths();
	virtual double findCorrelation(MtzFFTPtr one, MtzFFTPtr two);
	virtual void calculate();
private:
	std::string _csv;

	static std::vector<std::string> _names;
	static std::map<std::string, int> _ids;
	static std::map<std::string, std::vector<double> > _vectors;
	std::vector<double> _averageVec;
	ClusterList *_list;
};

#endif
