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
	CrystFELInput(std::string stream, std::string geom, 
	              std::string spg, double res, int max)
{
	_max = max;
	_stream = stream;
	_geom = geom;
	_spg = spg;
	_res = res;
}

	Group *process();
private:
	std::vector<struct image> loadStream();

	std::string _stream;
	std::string _geom;
	std::string _spg;
	double _res;
	int _max;

	struct detector *_detector;
};

#endif
