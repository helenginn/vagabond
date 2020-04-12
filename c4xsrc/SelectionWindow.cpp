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
#include "GLPoint.h"
#include "KeeperGL.h"

SelectionWindow::SelectionWindow(QWidget *parent, KeeperGL *keeper)
: QGraphicsView(parent)
{
	_keeper = keeper;
	setGeometry(0, 0, parent->width(), parent->height());
	QBrush brush(Qt::transparent);
	setBackgroundBrush(brush);
	setStyleSheet("background-color: transparent;");
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	_selectMode = false;
	_removeMode = false;
	_startX = 0;
	_startY = 0;
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
		return;
	}

	int endX = e->x();
	int endY = e->y();
	
	float x1 = _startX / (float)scene()->width();
	float y1 = endY / (float)scene()->height();
	float y2 = _startY / (float)scene()->height();
	float x2 = endX / (float)scene()->width();

	x1 = 2 * x1 - 1;
	x2 = 2 * x2 - 1;
	y1 = 2 * y1 - 1;
	y2 = 2 * y2 - 1;
	
	y1 *= -1;
	y2 *= -1;
	
	std::cout << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
	
	int add = 1;
	
	if (_startX == endX && _startY == endY)
	{
		std::cout << "Equivalent" << std::endl;
		add = 0;
	}
	
	if (_removeMode)
	{
		add = -1;
	}

	_points->selectInWindow(x1, y1, x2, y2, add);

	scene()->clear();
}

void SelectionWindow::mouseMoveEvent(QMouseEvent *e)
{
	if (_selectMode == false && _removeMode == false)
	{
		_keeper->mouseMoveEvent(e);
		return;
	}
	
	QRectF r;
	r.setTopLeft(QPoint(_startX, _startY));
	r.setBottomRight(QPoint(e->x(), e->y()));
	
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
	
	_startX = e->x();
	_startY = e->y();
}

void SelectionWindow::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Shift)
	{
		_selectMode = true;
	}
	if (event->key() == Qt::Key_Control)
	{
		_removeMode = true;
	}
}

void SelectionWindow::keyReleaseEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Shift)
	{
		_selectMode = false;
	}
	if (event->key() == Qt::Key_Control)
	{
		_removeMode = false;
	}
}

