// Fuck COV
// Copyright (C) 2017-2018 Helen Ginn
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

#include "Averager.h"
#include "KeeperGL.h"
#include "GLAxis.h"
#include "GLPoint.h"
#include "HKLView.h"
#include "MtzFFT.h"
#include "MtzFile.h"
#include <QApplication>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QWindow>
#include <QTimer>
#include <iostream>
#include <libsrc/mat4x4.h>

#define MOUSE_SENSITIVITY 500
#define PAN_SENSITIVITY 30
#define START_Z 0

void KeeperGL::initializeGL()
{
	_controlPressed = false;
	_lastX = 0;
	_lastY = 0;
	_moving = false;
	initializeOpenGLFunctions();

	glClearColor(1.0, 1.0, 1.0, 1.0);

//	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_BLEND);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	setupCamera();

	initialisePrograms();
}

void KeeperGL::updateProjection()
{
	zNear = 0.1;
	zFar = 50;

	double side = 0.5;
	float aspect = height() / width();
	_projMat = mat4x4_ortho(side, -side, side * aspect, -side * aspect,
	                       zNear, zFar);
}

void KeeperGL::updateCamera(void)
{
	vec3 alteration = empty_vec3();
	_totalCentroid = empty_vec3();
	
	vec3 centre = _centre;
	centre.x = 0;
	centre.y = 0;
	
	vec3 negCentre = centre;
	vec3_mult(&centre, -1);

	mat4x4 change = make_mat4x4();
	mat4x4_translate(&change, negCentre);
	mat4x4_rotate(&change, camAlpha, camBeta, camGamma);
	mat4x4_translate(&change, centre);

	mat4x4 transMat = make_mat4x4();
	_centre = vec3_add_vec3(_centre, _translation);
	_translation = vec3_add_vec3(_translation, alteration);
	mat4x4_translate(&transMat, _translation);

	mat4x4 tmp = mat4x4_mult_mat4x4(change, transMat);
	_model = mat4x4_mult_mat4x4(_model, tmp);
	
	mat3x3 rot = mat4x4_get_rot(_model);
	double det = mat3x3_determinant(rot);
	
	mat4x4_mult_scalar(&_model, _scale);
	_scale = 1;
	
	if (det < 0.05 && _store != NULL)
	{
		_model = make_mat4x4();
	}

	camAlpha = 0; camBeta = 0; camGamma = 0;
	_translation = make_vec3(0, 0, 0);
	
	if (_store != NULL)
	{
		*_store = _model;
		std::cout << "Store now" << std::endl;
		std::cout << mat4x4_desc(*_store) << std::endl;
	}
	
	for (int i = 0; i < 16; i++)
	{
		if (_model.vals[i] != _model.vals[i])
		{
			_model = make_mat4x4();
			std::cout << "Resetting model matrix" << std::endl;
		}
	}
}

void KeeperGL::setupCamera(void)
{
	_translation = make_vec3(0, 0, START_Z);
	_transOnly = make_vec3(0, 0, 0);
	_totalCentroid = make_vec3(0, 0, 0);
	_centre = make_vec3(0, 0, 0);
	camAlpha = 0;
	camBeta = 0;
	camGamma = 0;
	_scale = 1;

	updateProjection();

	updateCamera();
}

KeeperGL::KeeperGL(QWidget *p) : QOpenGLWidget(p)
{
	_store = NULL;
	_model = make_mat4x4();
	_rotMat = make_mat4x4();
	setGeometry(0, 0, p->width(), p->height());
	_points = NULL;
	_hklView = NULL;
}

std::vector<MtzFFTPtr> KeeperGL::getMtzs()
{
	return _ave->mtzs();
}

std::vector<MtzFFTPtr> KeeperGL::getMtzsFromSelection()
{
	std::vector<MtzFFTPtr> mtzs;
	for (size_t i = 0; i < _ave->mtzCount(); i++)
	{
		MtzFile *file = _ave->getMtz(i)->getMtzFile();
		if (file->isSelected())
		{
			mtzs.push_back(_ave->getMtz(i));
		}
	}
	
	return mtzs;
}

void KeeperGL::addAxes()
{
	GLAxis *axis = new GLAxis(make_vec3(1, 0, 0));
	_axes.push_back(axis);
	_renderMe.push_back(axis);
	axis = new GLAxis(make_vec3(0, 1, 0));
	_axes.push_back(axis);
	_renderMe.push_back(axis);
	axis = new GLAxis(make_vec3(0, 0, 1));
	_axes.push_back(axis);
	_renderMe.push_back(axis);
	
	
	update();
}

void KeeperGL::addSVDPoints(Averager *ave)
{
	_ave = ave;
	_points = new GLPoint();
	_points->setKeeper(this);
	_points->setAverager(ave);
	_renderMe.push_back(_points);
}

void KeeperGL::addHKLView(VagFFTPtr fft, double scale)
{
	_hklView = new HKLView(fft, scale);
	_renderMe.push_back(_hklView);
	_hklView->setKeeper(this);
	_hklView->repopulate();
	updateCamera();
	update();
}

void KeeperGL::preparePanels(int n)
{
	_axes.reserve(n);
}

void KeeperGL::panned(double x, double y)
{

}

void KeeperGL::draggedLeftMouse(double x, double y)
{
	x /= MOUSE_SENSITIVITY;
	y /= MOUSE_SENSITIVITY;

	camAlpha += y;
	camBeta += x;

	updateCamera();
	update();
}

void KeeperGL::draggedRightMouse(double x, double y)
{
	_scale *= (1 - y / MOUSE_SENSITIVITY);
	updateCamera();
	update();
}

void KeeperGL::mousePressEvent(QMouseEvent *e)
{
	_lastX = e->x();
	_lastY = e->y();
	_mouseButton = e->button();
	_moving = false;
}

void KeeperGL::mouseMoveEvent(QMouseEvent *e)
{
	if (_mouseButton == Qt::NoButton)
	{
		return;
	}

	_moving = true;

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
			panned(xDiff * 2, yDiff * 2);
		}
		else
		{
			draggedLeftMouse(-xDiff * 4, -yDiff * 4);
		}
	}
	else if (_mouseButton == Qt::RightButton)
	{
		draggedRightMouse(xDiff * PAN_SENSITIVITY, yDiff * PAN_SENSITIVITY);
	}
}
void KeeperGL::initialisePrograms()
{
	for (unsigned int i = 0; i < _renderMe.size(); i++)
	{
		_renderMe[i]->initialisePrograms();
	}
}

void KeeperGL::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	for (unsigned int i = 0; i < _renderMe.size(); i++)
	{
		_renderMe[i]->render(this);
	}
}

void KeeperGL::setAverager(Averager *ave)
{
	_ave = ave;
	_points->setAverager(ave);
}
