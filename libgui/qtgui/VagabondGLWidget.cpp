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
	_shiftPressed = false;
	setFocus();
}

void VagabondGLWidget::keyPressEvent(QKeyEvent *event)
{
	Density2GLPtr active = keeper->activeDensity();
	
	if (event->key() == Qt::Key_Alt || event->key() == Qt::Key_Control)
	{
		_controlPressed = true;
	}
	else if (event->key() == Qt::Key_Shift)
	{
		_shiftPressed = true;
		keeper->setAdding(true);
	}
	else if (event->key() == Qt::Key_D)
	{
		keeper->toggleVisibleDensity();
	}
	else if (active && (event->key() == Qt::Key_Plus ||
	                    event->key() == Qt::Key_Equal))
	{
		active->nudgeDensity(1);
	}
	else if (event->key() == Qt::Key_Minus && active)
	{
		active->nudgeDensity(-1);
	}
	else if (event->key() == Qt::Key_B)
	{
		keeper->toggleBondView();
	}
	else if (event->key() == Qt::Key_L)
	{
		_vag->toggleLog();
	}
	else if (event->key() == Qt::Key_R)
	{
		if (!keeper->isRefiningManually())
		{
			/* start it */
			_vag->setInstructionType(InstructionTypeManualRefine);
			_vag->wakeup();
		}
		else
		{
			/* send instruction to stop asap */
			keeper->cancelRefine();
		}
	}
	else if (event->key() == Qt::Key_Space)
	{
		keeper->focusOnSelected();
	}
	else if (event->key() == Qt::Key_S)
	{
		keeper->splitSelected();
	}
	else if (event->key() == Qt::Key_X)
	{
		keeper->deleteSelected();
	}
	else if (event->key() == Qt::Key_K)
	{
		keeper->toggleKicks();
	}
	else if (event->key() == Qt::Key_G)
	{
		_vag->gotoResidueDialogue();
	}
	else if (event->key() == Qt::Key_Comma)
	{
		keeper->advanceMonomer(-1);
	}
	else if (event->key() == Qt::Key_Period)
	{
		keeper->advanceMonomer(1);
	}
}

void VagabondGLWidget::keyReleaseEvent(QKeyEvent *event)
{
	_controlPressed = false;
	_shiftPressed = false;
	keeper->setAdding(false);
}

void VagabondGLWidget::mousePressEvent(QMouseEvent *e)
{
	_lastX = e->x();
	_lastY = e->y();
	_mouseButton = e->button();
	_moving = false;
	
	if (keeper->isRefiningManually() && e->button() == Qt::RightButton)
	{
		double x = e->x(); double y = e->y();
		convertCoords(&x, &y);
		keeper->setModelRay(x, y);
		keeper->setMouseRefine(true);
	}
}

void VagabondGLWidget::convertCoords(double *x, double *y)
{
	double w = width();
	double h = height();

	*x = 2 * *x / w - 1.0;
	*y =  - (2 * *y / h - 1.0);
}

void VagabondGLWidget::mouseReleaseEvent(QMouseEvent *e)
{
	if (keeper->isRefiningManually() && e->button() == Qt::RightButton)
	{
		keeper->setMouseRefine(false);
		return;
	}
	
	if (!_moving)
	{
		// this was just a click
		double prop_x = _lastX;
		double prop_y = _lastY;
		convertCoords(&prop_x, &prop_y);
		
		if (keeper)
		{
			keeper->findAtomAtXY(prop_x, prop_y);
		}
	}
	
	_mouseButton = Qt::NoButton;
}

void VagabondGLWidget::mouseMoveEvent(QMouseEvent *e)
{
	if (_mouseButton == Qt::NoButton)
	{
		return;
	}

	_moving = true;
	
	if (keeper->isRefiningManually() && e->button() == Qt::RightButton)
	{
		double prop_x = e->x();
		double prop_y = e->y();

		convertCoords(&prop_x, &prop_y);
		keeper->setModelRay(prop_x, prop_y);
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
			keeper->draggedLeftMouse(-xDiff * 4, -yDiff * 4);
		}
	}
	else if (_mouseButton == Qt::RightButton &&
	         !keeper->isRefiningManually())
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
	if (keeper)
	{
		keeper->getDensity2GL()->makeNewDensity(crystal);
		keeper->getWire2GL()->makeNewDensity(crystal);
		keeper->getDiffDens2GL()->makeNewDensity(crystal);
	}
}

VagabondGLWidget::~VagabondGLWidget()
{
	keeper->cleanup();
	delete keeper;
	delete timer;
}

void VagabondGLWidget::manualRefine()
{
	keeper->manualRefine();
}
