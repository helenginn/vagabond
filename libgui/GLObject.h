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
#include "Picture.h"
#include <vector>
#include "../libsrc/mat4x4.h"
#include <iostream>
#include <QOpenGLDebugLogger>

vec3 vec_from_pos(GLfloat *pos);

inline void checkErrors()
{
	GLenum err = 0;//glGetError();

	if (err != 0)
	{
		std::cout << "OUCH!" << std::endl;
	}
}

class GLObject : protected QOpenGLFunctions
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
	
	virtual void pause(bool on) {};

	vec3 fixCentroid(vec3 newCentroid);

	void setFocalPoint(vec3 vec);

	vec3 transformPosByModel(vec3 pos)
	{
		vec3 newPos = mat4x4_mult_vec(modelMat, pos);
		return newPos;
	}

	virtual bool isVagabond2GL()
	{
		return false;
	}
protected:
	mat4x4 projMat;
	mat4x4 modelMat;

	std::vector<Vertex> _vertices;
	std::vector<GLuint> _indices;

	void rebindProgram();
	virtual void bindTextures() {};
	
	void bindOneTexture(Picture &pic);
	
	void clearVertices()
	{
		_vertices.clear();
		_indices.clear();
	}
	
	bool _extra;
	bool _usesLighting;
	bool _usesFocalDepth;
	GLfloat _lightPos[3];
	GLfloat _xAxis[3];
	GLfloat _focalPos[3];

	bool _backToFront;
	GLenum _renderType;
	std::string *_fragShader;
	std::string *_vertShader;
	
	std::vector<GLuint> _textures;
	void reorderIndices();
private:
	void calculateCentroid();

	std::vector<IndexTrio> _temp; // stores with model mat
	static bool index_behind_index(IndexTrio one, IndexTrio two);
	static bool index_in_front_of_index(IndexTrio one, IndexTrio two);

	GLint _projectionUniform;
	GLint _modelUniform;
	GLint _lightUniform;
	GLint _focalUniform;
	GLint _xUniform;
	GLuint _colorRenderBuffer;
	GLuint _depthRenderBuffer;
	GLuint _vbo;
	GLuint _bufferID;

	GLuint _colorSlot;
	GLuint _program;
	vec3 _centroid;
	vec3 _focus;
};

#endif /* GLObject_hpp */
