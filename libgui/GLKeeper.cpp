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
#include "Bonds2GL.h"
#include "Atoms2GL.h"
#include "../libsrc/vec3.h"

bool GLKeeper::everMovedMouse = false;

void GLKeeper::updateProjection()
{
	zNear = 10;
	zFar = 100;

	float correctedNear = zNear;
	if (zNear <= 0.1) correctedNear = 0.1;
	double side = 0.5;
	float aspect = height / width;
	projMat = mat4x4_frustum(side, -side, side * aspect, -side * aspect,
	                         correctedNear, zFar);
}

void GLKeeper::setupCamera(void)
{
	_translation = make_vec3(0, 0, START_Z);
	_transOnly = make_vec3(0, 0, 0);
	_totalCentroid = make_vec3(0, 0, 0);
	_centre = make_vec3(0, 0, 0);
	camAlpha = 0;
	camBeta = 0;
	camGamma = 0;
	modelMat = make_mat4x4();
	rotMat = make_mat4x4();

	updateProjection();

	updateCamera();
}

void GLKeeper::updateCamera(void)
{
	vec3 alteration = _objects[0]->fixCentroid(_totalCentroid);
	_totalCentroid = _objects[0]->getCentroid();

	vec3 centre = _centre;
	centre.x = 0;
	centre.y = 0;
//	mat4x4 inv = mat4x4_inverse(modelMat);
//	centre = mat4x4_mult_vec(inv, centre);
	
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
	modelMat = mat4x4_mult_mat4x4(modelMat, tmp);

	camAlpha = 0; camBeta = 0; camGamma = 0;
	_translation = make_vec3(0, 0, 0);
}

void GLKeeper::focusOnPosition(vec3 pos)
{
	vec3 newPos = transformPosByModel(pos);
	_centre = vec3_add_vec3(_translation, newPos);
	vec3_mult(&newPos, -1);
	newPos.z -= 13;

	vec3 diff = vec3_subtract_vec3(newPos, _translation);
	_translation = vec3_add_vec3(_translation, newPos);
}

GLKeeper::GLKeeper(int newWidth, int newHeight)
{
	_setup.lock();

	width = newWidth;
	height = newHeight;
	_centre = empty_vec3();
	_densityState = 1;

	#ifdef SETUP_BUFFERS
	setupBuffers();
	#endif // SETUP_BUFFERS

	/* Bond model render */
	_allBond2GL = Bonds2GLPtr(new Bonds2GL(false));
	_totalCentroid = _allBond2GL->getCentroid();
	
	/* Average pos render */
	_aveBond2GL = Bonds2GLPtr(new Bonds2GL(true));
	_aveBond2GL->setEnabled(false);
	
	/* Atom pos render */
	_atoms2GL = Atoms2GLPtr(new Atoms2GL());

	/* Density render */
	_density2GL = Density2GLPtr(new Density2GL());
	_density2GL->setKeeper(this);
	_density2GL->recalculate();

	/* Difference density render */
	_diffDens2GL = Density2GLPtr(new Density2GL());
	_diffDens2GL->setKeeper(this);
	_diffDens2GL->setDiffDensity(true);
	_diffDens2GL->setVisible(false);
	_diffDens2GL->recalculate();

	_objects.push_back(_allBond2GL);
	_objects.push_back(_aveBond2GL);
	_objects.push_back(_atoms2GL);
	_objects.push_back(_density2GL);
	_objects.push_back(_diffDens2GL);

	setupCamera();

	initialisePrograms();
	
	_setup.unlock();
}

Density2GLPtr GLKeeper::activeDensity()
{
	if (_densityState == 0)
	{
		return Density2GLPtr();
	}
	else if (_densityState == 1)
	{
		return getDensity2GL();
	}
	
	return getDiffDens2GL();
}

void GLKeeper::toggleVisibleDensity()
{
	if (_densityState == 0)
	{
		_densityState++;
		getDensity2GL()->setVisible(true);
		getDiffDens2GL()->setVisible(false);
	}
	else if (_densityState == 1)
	{
		_densityState++;
		getDensity2GL()->setVisible(true);
		getDiffDens2GL()->setVisible(true);
	}
	else
	{
		_densityState = 0;	
		getDensity2GL()->setVisible(false);
		getDiffDens2GL()->setVisible(false);
	}
}

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
		_objects[i]->setProjMat(projMat);
		_objects[i]->setModelMat(modelMat);
		_objects[i]->render();
	}
	
	_setup.unlock();
}

void GLKeeper::toggleBondView()
{
	bool enabled = _allBond2GL->isEnabled();
	_aveBond2GL->setEnabled(enabled);
	_allBond2GL->setEnabled(!enabled);
}

void GLKeeper::cleanup(void)
{

}

void GLKeeper::keyPressed(char key)
{

}

void GLKeeper::zoom(float x, float y, float z)
{
	_translation.x += x;
	_translation.y += y;
	_translation.z += z;
	
	_transOnly.x += x;
	_transOnly.y += y;
	_transOnly.z += z;

}

void GLKeeper::rotateAngles(float alpha, float beta, float gamma)
{
	camAlpha -= alpha;
	camBeta += beta;
	camGamma -= gamma;
}

void GLKeeper::panned(float x, float y)
{
	zoom(x / MOUSE_SENSITIVITY * 4, y / MOUSE_SENSITIVITY * 4, 0);
}

void GLKeeper::draggedLeftMouse(float x, float y)
{
	if (!everMovedMouse)
	{
		everMovedMouse = true;
		return;
	}

	rotateAngles(y / MOUSE_SENSITIVITY, x / MOUSE_SENSITIVITY, 0);
}

void GLKeeper::draggedRightMouse(float x, float y)
{
	zoom(0, 0, -y / MOUSE_SENSITIVITY * 10);
}

void GLKeeper::changeSize(int newWidth, int newHeight)
{
	width = newWidth;
	height = newHeight;

	updateProjection();
}

void GLKeeper::initialisePrograms()
{
	for (int i = 0; i < _objects.size(); i++)
	{
		_objects[i]->initialisePrograms();
	}
}

