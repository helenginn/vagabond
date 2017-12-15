//
//  GLKeeper.cpp
//  RaddoseViewer
//
//  Created by Helen Ginn on 18/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#define MOUSE_SENSITIVITY 500

#define FIELD_OF_VIEW 20
#define DEGREES_TO_RADIANS M_PI / 180 *

#define BACKGROUND_COLOR_RED 0
#define BACKGROUND_COLOR_GREEN 0
#define BACKGROUND_COLOR_BLUE 0

//#define SETUP_BUFFERS

#include "GLKeeper.h"
#include "vec3.h"

bool GLKeeper::everMovedMouse = false;

void GLKeeper::setupCamera(void)
{
	_translation = make_vec3(0, 0, 0);
	camAlpha = 0;
	camBeta = 0;
	camGamma = 0;
	_centre = make_vec3(0, 0, START_Z);
    zNear = 3;
    zFar = 50;
	//	modelMat.vals[11] -= centreZ;
	modelMat = make_mat4x4();
	rotMat = make_mat4x4();

	float correctedNear = zNear;
	if (zNear <= 0.1) correctedNear = 0.1;
	double side = 0.3;
	float aspect = height / width;
	float fieldSize = correctedNear * tanf(DEGREES_TO_RADIANS(FIELD_OF_VIEW) / 2.0);
	projMat = mat4x4_frustum(side, -side, side * aspect, -side * aspect, correctedNear, zFar);

	updateCamera();
}

void GLKeeper::updateCamera(void)
{
	vec3 alteration = _objects[0]->fixCentroid(_centre);
	_translation = vec3_add_vec3(_translation, alteration);
	_centre = vec3_subtract_vec3(_centre, alteration);

	vec3 centre = make_vec3(0, 0, 5);
	vec3 negCentre = centre;
	vec3_mult(&negCentre, -1);

	mat4x4 change = make_mat4x4();
	mat4x4_translate(&change, centre);
	mat4x4_rotate(&change, camAlpha, camBeta, camGamma);
	mat4x4_translate(&change, negCentre);

	mat4x4 transMat = make_mat4x4();
	mat4x4_translate(&transMat, _translation);

	modelMat = mat4x4_mult_mat4x4(modelMat, transMat);
	modelMat = mat4x4_mult_mat4x4(change, modelMat);

	camAlpha = 0; camBeta = 0; camGamma = 0;
	_translation = make_vec3(0, 0, 0);
}


GLKeeper::GLKeeper(int newWidth, int newHeight)
{
    width = newWidth;
    height = newHeight;
	_rendered = false;

#ifdef SETUP_BUFFERS
	setupBuffers();
#endif // SETUP_BUFFERS

	Vagabond2GLPtr v2gl = Vagabond2GLPtr(new Vagabond2GL());
	GLObjectPtr object = boost::static_pointer_cast<GLObject>(v2gl);

	_objects.push_back(object);

	initialisePrograms();
	setupCamera();

	glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void GLKeeper::render(void)
{
	glClearColor(0.9, 0.9, 0.9, 1.0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	updateCamera();

//	glViewport(0, 0, width, height);

	for (int i = 0; i < _objects.size(); i++)
	{
		_objects[i]->setProjMat(projMat);
		_objects[i]->setModelMat(modelMat);
		_objects[i]->render();
	}
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
	setShouldRender();

	zoom(0, 0, -y / MOUSE_SENSITIVITY * 10);
}

void GLKeeper::changeSize(int newWidth, int newHeight)
{
    width = newWidth;
    height = newHeight;

	setShouldRender();
}

void GLKeeper::initialisePrograms()
{
	for (int i = 0; i < _objects.size(); i++)
	{
		_objects[i]->initialisePrograms();
	}
}

/*
void GLKeeper::setupBuffers(void)
{
    glGenRenderbuffers(1, &_depthRenderBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, _depthRenderBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, width, height);
    
    glGenRenderbuffers(1, &_colorRenderBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, _colorRenderBuffer);
    
    GLuint framebuffer;
    glGenFramebuffers(1, &framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, _colorRenderBuffer);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, _depthRenderBuffer);
}*/

