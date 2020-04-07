// Vagabond
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

/* holds a pointer and will provide getters and setters for any
 * double of your choice */

class Any
{
public:
	Any(double *ptr)
	{
		_ptr = ptr;
	}
	
	static double get(void *object)
	{
		return *(static_cast<Any *>(object)->_ptr);
	}
	
	static void set(void *object, double val)
	{
		*(static_cast<Any *>(object)->_ptr) = val;
	}

private:
	double *_ptr;

};

