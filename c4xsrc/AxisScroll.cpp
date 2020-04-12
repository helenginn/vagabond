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

#include "AxisScroll.h"
#include "Averager.h"
#include "GLPoint.h"
#include <QLabel>
#include <libsrc/FileReader.h>
#include <QPushButton>
#include <iostream>

AxisScroll::AxisScroll(QWidget *parent) : QScrollArea(parent)
{
	_points = NULL;
	_ave = NULL;

	setGeometry(0, 0, parent->width(), 60);
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
}

void AxisScroll::makeLayout()
{
	if (_ave == NULL)
	{
		std::cout << "!!!" << std::endl;
		return;
	}
	
	_viewport = new QWidget(this);
	_viewport->setGeometry(0, 0, 40 * (_ave->mtzCount() + 1), 60);
	
	QLabel *l = new QLabel("Axes:", _viewport);
	l->setGeometry(0, 0, 60, 40);
	l->show();

	int left = 60;
	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		QPushButton *b = new QPushButton(_viewport);
		double w = _ave->getDiagW(i);

		b->setText(QString::fromStdString(f_to_str(w, 1)));
		b->setCheckable(true);
		b->setGeometry(left, 0, 60, 40);
		b->show();
		left += 60;
		
		_axes.push_back(b);
		connect(b, &QPushButton::clicked, this, &AxisScroll::pressed);
		
		if (i < 3)
		{
			b->setChecked(true);
		}
	}
	
	setWidget(_viewport);
}

void AxisScroll::pressed()
{
	int total = 0;
	std::vector<int> checked;

	for (size_t i = 0; i < _axes.size(); i++)
	{
		QPushButton *b = _axes[i];
		
		if (b->isChecked())
		{
			total++;
			checked.push_back(i);
		}
	}
	
	if (total != 3)
	{
		return;
	}
	
	_points->setAxes(checked[0], checked[1], checked[2]);
	_points->repopulate();
}
