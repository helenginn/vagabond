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

#ifndef __cluster4x__crystfel__
#define __cluster4x__crystfel__

#include <string>
#include <vector>

class Group;
struct detector;

class CrystFELInput
{
public:
	CrystFELInput(std::string streams, std::string geom, 
	              std::string spg, double res);
	
	void setSkipAndMax(int skip, int max)
	{
		_max = max;
		_skip = skip;
	}
	
	void onlyLoad(std::string preload)
	{
		_preload = preload;
	}

	Group *process();
private:
	std::vector<struct image> loadStream(std::string str);

	std::vector<std::string> _streams;
	std::string _geom;
	std::string _preload;
	std::string _spg;
	double _res;
	int _max;
	int _skip;

	struct detector *_detector;
};

#endif
