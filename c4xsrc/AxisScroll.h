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

#ifndef __fuck__cov__AxisScroll__
#define __fuck__cov__AxisScroll__

#include <QScrollArea>

class Averager;
class GLPoint;
class QPushButton;
class QWidget;

class AxisScroll : public QScrollArea
{
Q_OBJECT
public:
	AxisScroll(QWidget *parent);

	void setAverager(Averager *ave)
	{
		_ave = ave;
	}

	void setPoints(GLPoint *points)
	{
		_points = points;
	}
	
	void makeLayout();
public slots:
	void pressed();

private:
	std::vector<QPushButton *> _axes;
	Averager *_ave;
	GLPoint *_points;
	QWidget *_viewport;

};


#endif
