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

#include "../subprojects/helen3d/libsrc/Frameworks.h"
#include "../subprojects/helen3d/libsrc/SlipGL.h"
#include "../libsrc/mat4x4.h"
#include <math.h>
#include <vector>
#include "Vagabond2GL.h"

class GLKeeper : public SlipGL
{
public:
	GLKeeper(QWidget *parent);
	void changeSize(int newWidth, int newHeight);

//	virtual void render(void);
//	void pause(bool on);

//	void rotateAngles(float alpha, float beta, float gamma);

//	void keyPressed(char key);
//	void draggedLeftMouse(float x, float y);
//	void draggedRightMouse(float x, float y);
//	void panned(float x, float y);
	
private:

};

#endif /* defined(__RaddoseViewer__GLKeeper__) */

