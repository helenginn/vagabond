//
//  Density2GL.h
//  VagabondViewer
//
//  Created by Helen Ginn on 22/7/2018.
//  Copyright © 2018 Ginn. All rights reserved.
//

#include <mutex>
#include "GLObject.h"
#include <vector>
#include "../libsrc/shared_ptrs.h"

struct ivec
{
	int x;
	int y;
	int z;
};

class GLKeeper;

class Density2GL : public GLObject
{
public:
	Density2GL()
	{
		_recalculate = 0;
		_renderType = GL_TRIANGLES;
		_resolution = 0.8;
		_cubeIndices = std::vector<std::vector<GLuint> >();
		_offset = make_vec3(-12, 3, -20);
		_visible = true;
		_threshold = 2.4;
		_backToFront = true;
	}
	
	virtual void render();
	
	void makeNewDensity(CrystalPtr crystal = CrystalPtr());
	void nudgeDensity(int dir);
	
	void toggleVisible()
	{
		_visible = !_visible;
	}
	
	void recalculate()
	{
		_recalculate = true;
	}
	
	void setKeeper(GLKeeper *keeper)
	{
		_keeper = keeper;
	}
private:
	CrystalPtr _crystal;
	GLKeeper *_keeper;
	void calculateContouring(CrystalPtr crystal);
	void makeUniformGrid();
	void setupIndexTable();
	int _recalculate;
	std::mutex _renderLock;
	
	vec3 _offset;
	vec3 getCentreOffset();
	double _resolution;
	double _threshold;
	ivec _dims;
	bool _visible;
	std::vector<std::vector<GLuint> > _cubeIndices;
};
