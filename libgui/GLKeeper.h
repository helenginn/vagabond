//
//  GLKeeper.h
//  RaddoseViewer
//
//  Created by Helen Ginn on 18/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#ifndef __RaddoseViewer__GLKeeper__
#define __RaddoseViewer__GLKeeper__

#include <iostream>

#include "Frameworks.h"
#include "../libsrc/mat4x4.h"
#include <math.h>
#include <vector>
#include "Vagabond2GL.h"

class GLKeeper : GLObject
{
private:
	std::vector<GLObjectPtr> _objects;

	float camAlpha, camBeta, camGamma;
	float zNear, zFar;
	static bool everMovedMouse;
    GLfloat width, height;
	bool _rendered;

	vec3 _centre;
	vec3 _translation;
	vec3 _totalCentroid;

    void setupVBOs (void);
    void setupBuffers(void);
    void setupCamera(void);
    void updateCamera(void);

	mat4x4 rotMat;

	void zoom(float x, float y, float z);

    GLubyte *newIndices;
	void initialisePrograms();

public:
    GLKeeper(int width, int height);
    void changeSize(int newWidth, int newHeight);

    virtual void render(void);
    void cleanup(void);
    
    void rotateAngles(float alpha, float beta, float gamma);
    
    void keyPressed(char key);
    void draggedLeftMouse(float x, float y);
	void draggedRightMouse(float x, float y);
	void panned(float x, float y);

	bool shouldRender()
	{
		return !_rendered;
	}

	void setShouldRender()
	{
		_rendered = false;
	}
};

#endif /* defined(__RaddoseViewer__GLKeeper__) */

