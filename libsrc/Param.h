// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
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

#ifndef __vagabond__param__
#define __vagabond__param__

#include <cstddef>

class Param
{
public:
	static void setValue(void *object, double value)
	{
		static_cast<Param *>(object)->_value = value;
	}
	
	static double getValue(void *object)
	{
		if (object == NULL)
		{
			 return 0;
		}
		return static_cast<Param *>(object)->_value;
	}
	
	double value()
	{
		return _value;
	}
	
	void set_value(double value)
	{
		_value = value;
	}

private:
	double _value;
};


#endif
