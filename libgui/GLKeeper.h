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

class GLKeeper : public GLObject
{
public:
	GLKeeper(int width, int height);
	void changeSize(int newWidth, int newHeight);
	void toggleBondView();
	void toggleVisibleDensity();
	Density2GLPtr activeDensity();

	void manualRefine();
	virtual void render(void);
	void cleanup(void);
	void focusOnPosition(vec3 pos, double dist = 13);   
	void pause(bool on);

	void rotateAngles(float alpha, float beta, float gamma);
	AtomPtr findAtomAtXY(double x, double y);

	void keyPressed(char key);
	void draggedLeftMouse(float x, float y);
	void draggedRightMouse(float x, float y);
	void panned(float x, float y);
	
	Vagabond2GLPtr getMulti2GL()
	{
		return _multi2GL;
	}
	
	Density2GLPtr getDiffDens2GL()
	{
		return _diffDens2GL;
	}
	
	Density2GLPtr getDensity2GL()
	{
		return _density2GL;
	}
	
	Density2GLPtr getWire2GL()
	{
		return _wire2GL;
	}
	
	Density2GLPtr getOrig2GL()
	{
		return _orig2GL;
	}
	
	Density2GLPtr getOrigDiff2GL()
	{
		return _origDiff2GL;
	}
	
	WarpGLPtr getWarpGL()
	{
		return _warpGL;
	}
	
	vec3 getCentre()
	{
		return _centre;
	}
	
	vec3 getTranslation()
	{
		return _transOnly;
	}
	
	bool isRefiningManually();
	void cancelRefine();
	void resetSelection();
	void setMouseRefine(bool val);
	
	void toggleKicks();
	void novalentSelected();
	void deleteSelected();
	void splitSelected();
	void focusOnSelected();
	void selectResidue(std::string chain, int number);
	void advanceMonomer(int dir);
	void setAdding(bool val);
	void setModelRay(double x, double y);
private:
	std::vector<GLObjectPtr> _objects;

	float camAlpha, camBeta, camGamma;
	float zNear, zFar;
	static bool everMovedMouse;
	GLfloat width, height;

	vec3 _centre;
	vec3 _translation;
	vec3 _transOnly;
	vec3 _totalCentroid;

	void setupVBOs (void);
	void setupBuffers(void);
	void setupCamera(void);
	void updateCamera(void);
	void updateProjection();
	void setFocalPoint(vec3 pos);
	void subtractFromFocalPoint(vec3 subtract);
    
	mat4x4 rotMat;

	void zoom(float x, float y, float z);

	GLubyte *newIndices;
	void initialisePrograms();
	GLObjectPtr activeObject();

	Density2GLPtr _density2GL;
	Density2GLPtr _wire2GL;
	Density2GLPtr _orig2GL;
	Density2GLPtr _origDiff2GL;
	Density2GLPtr _diffDens2GL;
	std::mutex _setup;

	int _densityState;
	
	WarpGLPtr _warpGL;
	Vagabond2GLPtr _allBond2GL;
	Vagabond2GLPtr _aveBond2GL;
	Vagabond2GLPtr _atoms2GL;
	Vagabond2GLPtr _multi2GL;
	Selected2GLPtr _selected2GL;
};

#endif /* defined(__RaddoseViewer__GLKeeper__) */

