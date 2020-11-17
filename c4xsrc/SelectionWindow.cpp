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

#include <QKeyEvent>
#include <iostream>
#include "SelectionWindow.h"
#include <Plot3D.h>
#include "KeeperGL.h"

SelectionWindow::SelectionWindow(QWidget *parent, KeeperGL *keeper)
: QGraphicsView(parent)
{
	_keeper = keeper;
	_plot = NULL;
	_selectMode = false;
	_removeMode = false;
	_startX = -1;
	_startY = -1;

	setGeometry(0, 0, parent->width(), parent->height());
	QBrush brush(Qt::transparent);
	setBackgroundBrush(brush);
	setStyleSheet("background-color: transparent;");
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	QGraphicsScene *scene = new QGraphicsScene(this);
	scene->setSceneRect(0, 0, width(), height());
	setScene(scene);
}

void SelectionWindow::resizeEvent(QResizeEvent *e)
{
	scene()->setSceneRect(0, 0, width(), height());
}

void SelectionWindow::mouseReleaseEvent(QMouseEvent *e)
{
	if (_selectMode == false && _removeMode == false)
	{
		_keeper->mouseReleaseEvent(e);
		return;
	}

	QPoint p = mapFromGlobal(e->globalPos());
	int endX = p.x();
	int endY = p.y();

	float x1 = std::min(_startX, endX);
	float y1 = std::max(_startY, endY);
	float x2 = std::max(_startX, endX);
	float y2 = std::min(_startY, endY);
	
	x1 /= (float)scene()->width();
	x2 /= (float)scene()->width();
	y1 /= (float)scene()->height();
	y2 /= (float)scene()->height();

	x1 = 2 * x1 - 1;
	x2 = 2 * x2 - 1;
	y1 = 2 * y1 - 1;
	y2 = 2 * y2 - 1;
	
	y1 *= -1;
	y2 *= -1;

	int add = 1;
	
	if (_startX == endX && _startY == endY)
	{
		add = 0;
	}
	
	if (_removeMode)
	{
		add = -1;
	}

	_plot->selectInWindow(x1, y1, x2, y2, add);

	scene()->clear();

	_startX = -1;
	_startY = -1;
	_selectMode = false;
	_removeMode = false;
}

void SelectionWindow::mouseMoveEvent(QMouseEvent *e)
{
	if (_selectMode == false && _removeMode == false)
	{
		_keeper->mouseMoveEvent(e);
		return;
	}
	
	QRectF r;

	QPoint p = mapFromGlobal(e->globalPos());
	
	if (_startX < 0 || _startY < 0)
	{
		return;
	}
	
	int left = std::min(_startX, p.x());
	int right = std::max(_startX, p.x());
	int top = std::min(_startY, p.y());
	int bottom = std::max(_startY, p.y());

	r.setTop(top);
	r.setBottom(bottom);
	r.setLeft(left);
	r.setRight(right);
	
	scene()->clear();
	scene()->addRect(r);
}

void SelectionWindow::mousePressEvent(QMouseEvent *e)
{
	if (_selectMode == false && _removeMode == false)
	{
		_keeper->mousePressEvent(e);
		return;
	}
	
	QPoint p = mapFromGlobal(e->globalPos());
	_startX = p.x();
	_startY = p.y();
}

void SelectionWindow::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Shift)
	{
		_selectMode = true;
	}
	else if (event->key() == Qt::Key_Control)
	{
		_removeMode = true;
	}
}

void SelectionWindow::keyReleaseEvent(QKeyEvent *event)
{
	/* don't change mode while dragging mouse */
	if (_startX >= 0 || _startY >= 0)
	{
		return;
	}

	if (event->key() == Qt::Key_Shift)
	{
		_selectMode = false;
	}
	if (event->key() == Qt::Key_Control)
	{
		_removeMode = false;
	}
}

