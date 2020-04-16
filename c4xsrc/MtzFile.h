// Clusterxxxx
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

#ifndef __clusterxxxx__MtzFile__
#define __clusterxxxx__MtzFile__

#include <string>

class MtzFile
{
public:
	MtzFile(std::string filename);

	std::string getFilename()
	{
		return _filename;
	}
	
	void setMarked(bool mark)
	{
		if (_dead) return;

		_mark = mark;
		_sele = false;
	}
	
	bool isMarked()
	{
		return _mark;
	}
	
	bool isDead()
	{
		return _dead;
	}
	
	void setDead(bool dead)
	{
		_dead = dead;
		_mark = false;
		_sele = false;
	}
	
	bool isSelected()
	{
		return _sele;
	}
	
	void flipSelected()
	{
		_sele = !_sele;
	}

	void setSelected(bool sele)
	{
		if (!_mark && !_dead)
		{
			_sele = sele;
		}
	}
private:
	std::string _filename;

	bool _mark;
	bool _sele;
	bool _dead;
};

#endif
