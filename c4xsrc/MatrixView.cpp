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

#include "MatrixView.h"
#include "MtzFFT.h"
#include "MtzFile.h"
#include "Group.h"
#include <QPainter>
#include <iostream>

MatrixView::MatrixView(Group *ave, int w, int h) : QImage(w + 1, h, QImage::Format_RGB32)
{
	_ave = ave;
	_contrast = 1;
	fill(Qt::white);
}

void MatrixView::populate()
{
	int num = _ave->mtzCount();
	double **raw = _ave->getRawPtr();

	QPainter painter(this);

	double box_size = ((double)width() / (double)(num));
	
	int red = 255;
	int green = 0;
	int blue = 0;

	for (int j = 0; j < num; j++)
	{
		MtzFile *file = _ave->getMtz(j)->getMtzFile();
		for (int i = 0; i < num; i++)
		{
			double val = raw[i][j];
			val /= _contrast;
			
			if (val > 2) val = 2;

			if (val != val) /* we go grey */
			{
				red = 100;
				green = 100;
				blue = 100;
			}
			else if (val <= 0) /* we go black */
			{
				val = std::min(-val, 1.);
				red = 0;
				green = 0;
				blue = 255 - val * 255;
			}
			else if (val < 0.5)
			{
				/* we go blue. */
				val = (0.5 - val ) * 2.;
				red = 255 - val * 255;
				green = 255 - val * 255;
				blue = 255;
			}
			else if (val >= 1.0) /* We go yellow. */
			{
				val -= 1; 
				red = 255;
				green = val * 255;
				blue = 0;
			}
			else if (val >= 0.5) /* We go red. */
			{
				val = (val - 0.5) * 2.0;
				red = 255;
				green = 255 - val * 255;
				blue = 255 - val * 255;
			}
			
			
			if (file->isDead())
			{
				red = 0;
				green = 0;
				blue = 0;
			}
			else if (file->isMarked())
			{
				red *= 0.5;
				green *= 0.5;
				blue *= 0.5;
			}
			else if (file->isSelected())
			{
				red = 255 - (255 - red) / 2;
				green = 255 - (255 - green) / 2;
				blue = 0;
			}

			QColor c = QColor(red, green, blue, 255);
			QPen p = QPen(c);
			QBrush b = QBrush(c, Qt::SolidPattern);
			painter.setPen(p);
			painter.setBrush(b);
			
			painter.drawRect(box_size * i, box_size * j,
			                 box_size + 1, box_size + 1);
		}
	}
}

void MatrixView::updateSelection()
{
	populate();
}
