// Cluster4x
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

#include <QMouseEvent>
#include <iostream>
#include "CorrelLabel.h"
#include "Screen.h"
#include "MtzFile.h"
#include "MatrixView.h"
#include "Group.h"
#include "ClusterList.h"

CorrelLabel::CorrelLabel(QWidget *parent, MatrixView *image,
Screen *screen) : QLabel(parent)
{
	_image = image;
	_screen = screen;
	_list = _screen->getList();
	_shift = false;
	_ctrl = false;
	_lastFile = -1;
	_startFile = -1;
}

void CorrelLabel::keyPressEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_Shift)
	{
		_shift = true;
	}
	if (e->key() == Qt::Key_Control)
	{
		_ctrl = true;
	}
}

void CorrelLabel::keyReleaseEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_Shift)
	{
		_shift = false;
	}
	if (e->key() == Qt::Key_Control)
	{
		_ctrl = false;
	}
	
	if (e->key() == Qt::Key_M)
	{
		_list->cycleCSV(true);
	}
	
	if (e->key() == Qt::Key_N)
	{
		_list->cycleCSV(false);
	}
}

void CorrelLabel::mouseReleaseEvent(QMouseEvent *e)
{

}

void CorrelLabel::mouseMoveEvent(QMouseEvent *e)
{
	if (!_shift && !_ctrl)
	{
		return;
	}

	QPoint p = mapFromGlobal(e->globalPos());
	int file = findMtzForXY(p.x(), p.y());
	if (file < 0)
	{
		return;
	}

	int dir = (file > _startFile) ? 1 : -1;
	
	Group *ave = _image->getGroup();
	bool sele = (_ctrl ? false : _shift);

	if (dir > 0 && file > _lastFile)
	{
		for (int i = _lastFile + 1; i <= file; i++)
		{
			ave->setMtzSelection(i, sele);
		}
	}
	else if (dir > 0 && file < _lastFile)
	{
		for (int i = _lastFile; i >= file; i--)
		{
			ave->setMtzSelection(i, !sele);
		}
	}
	else if (dir < 0 && file > _lastFile)
	{
		for (int i = _lastFile; i < file; i++)
		{
			ave->setMtzSelection(i, !sele);
		}
	}
	else if (dir < 0 && file < _lastFile)
	{
		for (int i = _lastFile - 1; i >= file; i--)
		{
			ave->setMtzSelection(i, sele);
		}
	}
	
	_lastFile = file;
	_screen->refreshSelection();
}

void CorrelLabel::mousePressEvent(QMouseEvent *e)
{
	if (!_shift && !_ctrl)
	{
		return;
	}

	QPoint p = mapFromGlobal(e->globalPos());
	int f = findMtzForXY(p.x(), p.y());
	
	if (f < 0)
	{
		return;
	}

	_lastFile = f;
	_startFile = f;

	bool sele = (_ctrl ? false : _shift);

	Group *ave = _image->getGroup();
	ave->setMtzSelection(f, sele);

	_screen->refreshSelection();
}

int CorrelLabel::findMtzForXY(int x, int y)
{
	Group *ave = _image->getGroup();
	int mtzs = ave->mtzCount();

	float ty = mtzs * y / (float)height();

	if (ty > mtzs)
	{
		return -1;
	}

	return ty;

}
