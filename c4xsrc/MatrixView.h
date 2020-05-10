// Fuck COV
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

#ifndef __fuck_cov__MatrixView__
#define __fuck_cov__MatrixView__

#include <QImage>
#include <map>

class Group;

class MatrixView : public QImage
{
public:
	MatrixView(Group *ave, int width = 500, int height = 500);

	void populate();
	
	void setContrast(double cont)
	{
		_contrast = cont;
	}
	
	void updateSelection();
	
	Group *getGroup()
	{
		return _ave;
	}
private:
	double _contrast;
	Group *_ave;

};

#endif
