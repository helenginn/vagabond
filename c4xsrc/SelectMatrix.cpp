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

#include "SelectMatrix.h"
#include <QLineEdit>
#include <QPushButton>
#include <iostream>
#include <libsrc/mat3x3.h>
#include "ClusterList.h"

SelectMatrix::SelectMatrix(QWidget *widget) : QMainWindow(widget)
{
	setGeometry(300, 300, 320, 320);
	setWindowTitle("Reindex + translate");
	const int size = 50;
	const int margin = 10;

	int top = size;
	int left = size;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int holder = 0;
			QLineEdit *line = new QLineEdit(this);
			line->setGeometry(left, top, size, size);
			line->setPlaceholderText(QString::number(holder));
			line->show();
			_entries.push_back(line);
			left += size + margin;
		}
		
		top += size + margin;
		left = size;
	}
	
	top = size;
	left = size + (size + margin) * 3;

	for (int i = 0; i < 3; i++)
	{
		int holder = 0;
		QLineEdit *line = new QLineEdit(this);
		line->setGeometry(left, top, size, size);
		line->setPlaceholderText(QString::number(holder));
		line->show();
		_trans.push_back(line);
		top += size + margin;
	}

	QPushButton *p = new QPushButton("Load", this);
	p->setGeometry(left, top, 80, 40);
	p->show();

	connect(p, &QPushButton::clicked, this, &SelectMatrix::load);
}

void SelectMatrix::load()
{
	mat3x3 reindex = make_mat3x3();
	vec3 vec = empty_vec3();

	for (size_t i = 0; i < _entries.size() && i < 9; i++)
	{
		QString val = _entries[i]->text();
		int num = val.toInt();
		reindex.vals[i] = num;
	}

	for (size_t i = 0; i < _trans.size() && i < 3; i++)
	{
		QString val = _trans[i]->text();
		double num = val.toDouble();
		*(&vec.x + i) = num;
	}
	
	std::cout << mat3x3_desc(reindex) << std::endl;
	std::cout << vec3_desc(vec) << std::endl;
	
	if (_list != NULL)
	{
		hide();
		deleteLater();
		_list->setReindexMatrix(reindex, vec);
	}
}
