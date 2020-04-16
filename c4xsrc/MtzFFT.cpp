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

#include "MtzFFT.h"
#include "MtzFile.h"

MtzFFT::MtzFFT(QTreeWidgetItem *parent, VagFFT &vag) : VagFFT(vag, -1), QTreeWidgetItem(parent)
{
	
}

MtzFFT::MtzFFT(QTreeWidgetItem *parent, MtzFFT &vag) : VagFFT(vag, -1), QTreeWidgetItem(parent)
{
	_file = vag._file;
	setText(0, vag.text(0));
}

void MtzFFT::updateText()
{
	QColor c = QColor(255, 255, 255, 255);
	if (getMtzFile()->isMarked())
	{
		c = QColor(255, 50, 50, 255);
	}
	else if (getMtzFile()->isSelected())
	{
		c = QColor(200, 200, 0, 255);
	}
	else if (getMtzFile()->isDead())
	{
		c = QColor(100, 100, 100, 255);
	}

	setBackground(0, QBrush(c));
}
