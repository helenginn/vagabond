//
//  GLKeeper.cpp
//  RaddoseViewer
//
//  Created by Helen Ginn on 18/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#define MOUSE_SENSITIVITY 500

#define DEGREES_TO_RADIANS M_PI / 180 *

#define BACKGROUND_COLOR_RED 0
#define BACKGROUND_COLOR_GREEN 0
#define BACKGROUND_COLOR_BLUE 0

#include "GLKeeper.h"
#include "Density2GL.h"
#include "Connect2GL.h"
#include "Multi2GL.h"
#include "Bonds2GL.h"
#include "Selected2GL.h"
#include <float.h>
#include "../libsrc/vec3.h"

/*
void GLKeeper::updateProjection()
{
	zNear = 5;
	zFar = 400;

	double side = 0.5;
	float aspect = height / width;
	_proj = mat4x4_frustum(-side, side, side * aspect, -side * aspect,
	                       zNear, zFar);
	_unproj = mat4x4_unfrustum(-side, side, side * aspect, -side * aspect,
	                       zNear, zFar);
}

void GLKeeper::setupCamera(void)
{
	_translation = make_vec3(0, 0, START_Z);
	_centre = make_vec3(0, 0, 0);
	camAlpha = 0;
	camBeta = 0;
	camGamma = 0;
	_model = make_mat4x4();
	rotMat = make_mat4x4();

	updateProjection();

	updateCamera();
}
*/

/*
void GLKeeper::updateCamera(void)
{
	vec3 centre = _centre;
	vec3 negCentre = _centre;
	vec3_mult(&negCentre, -1);

	mat4x4 change = make_mat4x4();
	mat4x4_translate(&change, negCentre);
	mat4x4_rotate(&change, camAlpha, camBeta, camGamma);
	mat4x4_translate(&change, centre);

	mat4x4 transMat = make_mat4x4();
	_centre = vec3_add_vec3(_centre, _translation);
	mat4x4_translate(&transMat, _translation);

	mat4x4 tmp = mat4x4_mult_mat4x4(change, transMat);
	_model = mat4x4_mult_mat4x4(tmp, _model);

	camAlpha = 0; camBeta = 0; camGamma = 0;
	_translation = make_vec3(0, 0, 0);
	
	_centre.x = 0;
	_centre.y = 0;
	setFocalPoint(negCentre);
}
*/

GLKeeper::GLKeeper(QWidget *parent) : SlipGL(parent)
{
}

/*
void GLKeeper::pause(bool on)
{
	for (int i = 0; i < _objects.size(); i++)
	{
		_objects[i]->pause(on);
	}
}

void GLKeeper::render(void)
{
	if (!_setup.try_lock())
	{
		return;
	}

	glClearColor(1.0, 1.0, 1.0, 1.0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	updateCamera();

	for (int i = 0; i < _objects.size(); i++)
	{
		_objects[i]->setProjMat(_proj);
		_objects[i]->setUnprojMat(_unproj);
		_objects[i]->setModelMat(_model);
		_objects[i]->render();
	}
	
	_setup.unlock();
}
*/

/*
void GLKeeper::keyPressed(char key)
{

}

void GLKeeper::zoom(float x, float y, float z)
{
	_translation.x += x;
	_translation.y += y;
	_translation.z += z;
}

void GLKeeper::rotateAngles(float alpha, float beta, float gamma)
{
	camAlpha -= alpha;
	camBeta += beta;
	camGamma -= gamma;
}

void GLKeeper::panned(float x, float y)
{
	zoom(-x / MOUSE_SENSITIVITY * 4, y / MOUSE_SENSITIVITY * 4, 0);
}

void GLKeeper::draggedLeftMouse(float x, float y)
{
	if (!everMovedMouse)
	{
		everMovedMouse = true;
		return;
	}

	rotateAngles(y / MOUSE_SENSITIVITY, -x / MOUSE_SENSITIVITY, 0);
}

void GLKeeper::draggedRightMouse(float x, float y)
{
	zoom(0, 0, -y / MOUSE_SENSITIVITY * 10);
}

void GLKeeper::initialisePrograms()
{
	for (int i = 0; i < _objects.size(); i++)
	{
		_objects[i]->initialisePrograms();
	}
}
*/

