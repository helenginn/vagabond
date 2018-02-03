//
//  GLObject.hpp
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef GLObject_hpp
#define GLObject_hpp

#include <stdio.h>
#include "Frameworks.h"
#include <vector>
#include "../libsrc/mat4x4.h"
#include <iostream>

const Vertex Vertices[] = {
	{{0.5, 0.0, -8}, {1, 0, 0}, {1, 0, 0, 1}},
	{{0.0, 0.5, -8}, {0, 1, 0}, {0.33, 0.33, 0.33, 1}},
	{{-0.5, 0, -8}, {0, 0, 1}, {0, 1, 0, 1}},
	{{0.0, -0.5, -8}, {0, 0, 1}, {0, 0, 1, 1}}
};

const GLuint Indices[] = {
	0, 1,
	1, 2,
	2, 3,
	3, 0
};

inline void checkErrors()
{
	GLenum err = glGetError();

	if (err != 0)
	{
		std::cout << "OUCH!" << std::endl;
	}
}

class GLObject
{
public:
	GLObject();
	virtual void render();
	void initialisePrograms(void);

	Vertex *vPointer()
	{
		return &_vertices[0];
	}

	size_t vSize()
	{
		return sizeof(Vertex) * _vertices.size();
	}

	GLuint *iPointer()
	{
		return &_indices[0];
	}

	size_t iSize()
	{
		return sizeof(GLuint) * _indices.size();
	}

	size_t indexCount()
	{
		return _indices.size();
	}

	void setProjMat(mat4x4 proj)
	{
		projMat = proj;
	}

	void setModelMat(mat4x4 model)
	{
		modelMat = model;
	}

	vec3 getCentroid()
	{
		calculateCentroid();
		return _centroid;
	}

	vec3 fixCentroid(vec3 newCentroid);

        vec3 transformPosByModel(vec3 pos)
        {
            vec3 newPos = mat4x4_mult_vec(modelMat, pos);
            return newPos;
        }
protected:
	mat4x4 projMat;
	mat4x4 modelMat;

	std::vector<Vertex> _vertices;
	std::vector<GLuint> _indices;

	void rebindProgram();
private:
	void calculateCentroid();

	GLint _projectionUniform;
	GLint _modelUniform;
	GLuint _colorRenderBuffer;
	GLuint _depthRenderBuffer;
	GLuint _vbo;
	GLuint _bufferID;

	GLuint _colorSlot;
	GLuint _program;
	vec3 _centroid;
};

#endif /* GLObject_hpp */
