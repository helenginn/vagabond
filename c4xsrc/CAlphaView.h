// 
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

#ifndef __cluster__calphaview__
#define __cluster__calphaview__

#include <libsrc/vec3.h>
#include "GLObject.h"

class Averager;
class MtzFile;

class CAlphaView : public QObject, public GLObject
{
Q_OBJECT
public:
	CAlphaView(MtzFile *mtz, vec3 centre = empty_vec3());
	CAlphaView(Averager *ave, vec3 centre = empty_vec3());

	void setKeeper(KeeperGL *gl)
	{
		_keeper = gl;
	}

	virtual void initialisePrograms();
	void repopulate();
	void recolour();
private:
	void addCAlpha(vec3 point);

	std::map<MtzFile *, size_t> _starts;
	std::map<MtzFile *, size_t> _ends;
	std::vector<MtzFile *> _mtzs;
	vec3 _centre;
	KeeperGL *_keeper;

};

#endif
