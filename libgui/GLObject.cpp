//
//  GLObject.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "GLObject.h"
#include "shader.h"
#include <iostream>

GLObject::GLObject()
{
	_projectionUniform = 0;
}

void GLObject::rebindProgram()
{
	glBufferData(GL_ARRAY_BUFFER, vSize(), vPointer(), GL_STATIC_DRAW);
	checkErrors();

	glBufferData(GL_ELEMENT_ARRAY_BUFFER, iSize(), iPointer(), GL_STATIC_DRAW);
	checkErrors();

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)(0 * sizeof(float)));
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)(3 * sizeof(float)));
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)(6 * sizeof(float)));
	checkErrors();

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
}

void GLObject::render()
{
	glUseProgram(_program);

	const char *uniform_name = "projection";
	_projectionUniform = glGetUniformLocation(_program, uniform_name);
	glUniformMatrix4fv(_projectionUniform, 1, GL_FALSE, &projMat.vals[0]);
	checkErrors();

	const char *model_name = "model";
	_modelUniform = glGetUniformLocation(_program, model_name);
	glUniformMatrix4fv(_modelUniform, 1, GL_FALSE, &modelMat.vals[0]);
	checkErrors();

	glDrawElements(GL_LINES, indexCount(), GL_UNSIGNED_INT, 0);

	glUseProgram(0);
}

void GLObject::initialisePrograms()
{
	GLint result;

	/* create program object and attach shaders */
	_program = glCreateProgram();
	shaderAttachFromFile(_program, GL_VERTEX_SHADER, "Shader.vsh");
	shaderAttachFromFile(_program, GL_FRAGMENT_SHADER, "Shader.fsh");

	glBindAttribLocation(_program, 0, "position");
	glBindAttribLocation(_program, 1, "normal");
	glBindAttribLocation(_program, 2, "color");
	glBindAttribLocation(_program, 3, "projection");
	_projectionUniform = 3;
	checkErrors();

	/* link the program and make sure that there were no errors */
	glLinkProgram(_program);
	glGetProgramiv(_program, GL_LINK_STATUS, &result);
	if (result == GL_FALSE)
	{
		std::cout << "sceneInit(): Program linking failed." << std::endl;

		/* delete the program */
		glDeleteProgram(_program);
		_program = 0;
	}
	glGenBuffers(1, &_bufferID);
	glGenBuffers(1, &_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, _bufferID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo);
	
	rebindProgram();
}

void GLObject::calculateCentroid()
{
	vec3 sum;

	for (int i = 0; i < _vertices.size(); i++)
	{
		vec3 add = make_vec3(_vertices[i].pos[0],
							 _vertices[i].pos[1],
							 _vertices[i].pos[2]);

		sum = vec3_add_vec3(sum, add);
	}

	double mult = 1 / (double)_vertices.size();
	vec3_mult(&sum, mult);

	if (!_vertices.size())
	{
		sum = make_vec3(0, 0, 0);
	}
	
	_centroid = sum;
}

vec3 GLObject::fixCentroid(vec3 newCentroid)
{
	if (!_vertices.size())
	{
		return make_vec3(0, 0, 0);
	}

	vec3 oldCentroid = getCentroid();
	vec3 addition = vec3_subtract_vec3(newCentroid, oldCentroid);

	return addition;
}
