// Fuck COV
// Copyright (C) 2020 Helen Ginn
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

#ifndef __fuck_cov__GLPoint__
#define __fuck_cov__GLPoint__

#include "GLObject.h"
#include "vec3.h"

class Group;
class KeeperGL;

class GLPoint : public QObject, public GLObject
{
Q_OBJECT
public:
	GLPoint();
	
	void setGroup(Group *ave);
	void addPoint(vec3 point);
	virtual void initialisePrograms();
	
	void setAxes(int a, int b, int c)
	{
		_a = a;
		_b = b;
		_c = c;
	}
	
	void setKeeper(KeeperGL *gl)
	{
		_keeper = gl;
	}

	/* add 1 = add to selection
	 * add 0 = replace selection
	 * add -1 = remove from selection */
	void selectInWindow(float x1, float y1, float x2, float y2,
	                    int add);
	void repopulate();
signals:
	void updateSelection();
private:
	Group *_ave;
	KeeperGL *_keeper;
	void recolour();
	
	int _a;
	int _b;
	int _c;
};

#endif

