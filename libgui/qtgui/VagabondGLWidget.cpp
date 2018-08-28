//
//  VagabondGLWidget.h
//  Vagabond
//
//  Created by Helen Ginn on 20/01/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include "VagabondGLWidget.h"
#include "../Density2GL.h"
#include "../../libsrc/shared_ptrs.h"

#define PAN_SENSITIVITY 30

VagabondGLWidget::VagabondGLWidget(QWidget *obj) : QOpenGLWidget(obj)
{
	keeper = NULL;
	timer = NULL;
	_mouseButton = Qt::NoButton;
	_lastX = 0; _lastY = 0;
	_controlPressed = false;
	setFocus();
}

void VagabondGLWidget::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Alt)
	{
		_controlPressed = true;
	}
	else if (event->key() == Qt::Key_D)
	{
		keeper->getDensity2GL()->toggleVisible();
	}
	else if (event->key() == Qt::Key_Plus ||
	         event->key() == Qt::Key_Equal)
	{
		keeper->getDensity2GL()->nudgeDensity(1);
	}
	else if (event->key() == Qt::Key_Minus)
	{
		keeper->getDensity2GL()->nudgeDensity(-1);
	}
}

void VagabondGLWidget::keyReleaseEvent(QKeyEvent *event)
{
	_controlPressed = false;
}

void VagabondGLWidget::mousePressEvent(QMouseEvent *e)
{
	_lastX = e->x();
	_lastY = e->y();
	_mouseButton = e->button();
}

void VagabondGLWidget::mouseReleaseEvent(QMouseEvent *e)
{
	_mouseButton = Qt::NoButton;
}

void VagabondGLWidget::mouseMoveEvent(QMouseEvent *e)
{
	if (_mouseButton == Qt::NoButton)
	{
		return;
	}

	double newX = e->x();
	double xDiff = _lastX - newX;
	double newY = e->y();
	double yDiff = _lastY - newY;
	_lastX = newX;
	_lastY = newY;

	if (_mouseButton == Qt::LeftButton)
	{
		if (_controlPressed)
		{
			keeper->panned(xDiff * 2, yDiff * 2);
		}
		else
		{
			keeper->draggedLeftMouse(xDiff * 4, yDiff * 4);
		}
	}
	else if (_mouseButton == Qt::RightButton)
	{
		keeper->draggedRightMouse(xDiff * PAN_SENSITIVITY,
		                          yDiff * PAN_SENSITIVITY);
	}
}

void VagabondGLWidget::initializeGL()
{
	keeper = new GLKeeper(width(), height());

	timer = new QTimer();
	timer->setInterval(30);
	connect(timer, SIGNAL(timeout()), this, SLOT(update()));
	timer->start();
}

void VagabondGLWidget::resizeGL()
{
	int w = width();
	int h = height();
	if (keeper)
	{
		keeper->changeSize(w, h);
	}
}

void VagabondGLWidget::paintGL()
{
	keeper->render();
}

void VagabondGLWidget::renderDensity(CrystalPtr crystal)
{
	keeper->getDensity2GL()->makeNewDensity(crystal);
}

VagabondGLWidget::~VagabondGLWidget()
{
	keeper->cleanup();
	delete keeper;
	delete timer;
}
