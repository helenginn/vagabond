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
#include <hcsrc/maths.h>
#include <QPainter>
#include <iostream>

MatrixView::MatrixView(Group *ave, int w, int h) : QImage(w, h, QImage::Format_RGB32)
{
	_ave = ave;
	_names = NULL;
	_contrast = 1;
	fill(Qt::white);
}

void MatrixView::populate()
{
	populateFromGroup();
}

std::string MatrixView::getName(int x, int y)
{
	if (x < 0 || y < 0 || x >= width() || y >= height())
	{
		return "";
	}

	double box_width = ((double)width() / (double)(_w));
	double box_height = ((double)height() / (double)(_h));

	int nx = x / box_width;
	int ny = y / box_height;

	return std::string(_names[nx][ny]);
}

void MatrixView::populateFromGroup()
{
	int num = _ave->mtzCount();
	double **raw = _ave->getRawPtr();
	QPainter painter(this);
	double box_size = ((double)width() / (double)(num));
	
	int red = 255;
	int green = 0;
	int blue = 0;

	for (int i = 0; i < num; i++)
	{
		MtzFile *file = _ave->getMtz(i)->getMtzFile();
		for (int j = 0; j < num; j++)
		{
			double val = raw[i][j];
			val /= _contrast;
			
			val_to_cluster4x_colour(val, &red, &green, &blue);
			
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
			
			painter.drawRect(box_size * j, box_size * i,
			                 box_size + 1, box_size + 1);
		}
	}
}

void MatrixView::populate(int w, int h, double **raw)
{
	_w = w;
	_h = h;

	QPainter painter(this);

	double box_width = ((double)width() / (double)(w));
	double box_height = ((double)height() / (double)(h));
	
	int red = 255;
	int green = 0;
	int blue = 0;

	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			double val = raw[i][j];
			val /= _contrast;
			
			val_to_cluster4x_colour(val, &red, &green, &blue);
			
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

			QColor c = QColor(red, green, blue, 255);
			QPen p = QPen(c);
			QBrush b = QBrush(c, Qt::SolidPattern);
			painter.setPen(p);
			painter.setBrush(b);
			
			painter.drawRect(box_width * i, box_height * j,
			                 box_width + 1, box_height + 1);
		}
	}
}

void MatrixView::updateSelection()
{
	populate();
}
