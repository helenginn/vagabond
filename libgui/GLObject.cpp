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
#include "Shaders/Shader_fsh.h"
#include "Shaders/Shader_vsh.h"

GLObject::GLObject()
{
    _projectionUniform = 0;
    initializeOpenGLFunctions();
	_renderType = GL_LINES;
	_fragShader = &Shader_fsh;
	_vertShader = &Shader_vsh;
	_extra = false;
	_backToFront = true;
	_usesLighting = false;
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

	if (_extra)
	{
		glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex),
		                      (void *)(10 * sizeof(float)));
	}

	if (_textures.size())
	{
		glVertexAttribPointer(4, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(14 * sizeof(float)));
	}

    checkErrors();

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

	if (_extra)
	{
		glEnableVertexAttribArray(3);
	}

	if (_textures.size())
	{
		glEnableVertexAttribArray(4); 
	}
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

	if (_usesLighting)
	{
		const char *light_name = "light_pos";
		_lightUniform = glGetUniformLocation(_program, light_name);
		glUniform3fv(_lightUniform, 1, &_lightPos[0]);
	}

	if (_textures.size())
	{
		glBindTexture(GL_TEXTURE_2D, _textures[0]);
	}
	
	if (_renderType == GL_POINTS)
	{
		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_POINT_SPRITE);
	}

    glDrawElements(_renderType, indexCount(), GL_UNSIGNED_INT, 0);

    glUseProgram(0);
}

void GLObject::initialisePrograms()
{
	GLint result;

	/* create program object and attach shaders */
	_program = glCreateProgram();

	Shader::shaderAttachFromFile(_program,  GL_FRAGMENT_SHADER, _fragShader->c_str(), true);
	Shader::shaderAttachFromFile(_program,  GL_VERTEX_SHADER, _vertShader->c_str(), true);

	glBindAttribLocation(_program, 0, "position");
	glBindAttribLocation(_program, 1, "normal");
	glBindAttribLocation(_program, 2, "color");

	if (!_extra)
	{
		glBindAttribLocation(_program, 3, "projection");
	}
	else
	{
		glBindAttribLocation(_program, 3, "extra");
	}
	
	if (_textures.size())
	{
		glBindAttribLocation(_program, 4, "tex");
	}

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
    
	bindTextures();
	
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

void GLObject::bindOneTexture(Picture &pic)
{
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pic.width, pic.height, 
	             0, GL_RGBA, GL_UNSIGNED_BYTE, pic.data);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	glGenerateMipmap(GL_TEXTURE_2D);
	checkErrors();
}

bool GLObject::index_behind_index(IndexTrio &one, IndexTrio &two)
{
	return (one.z > two.z);
}

bool GLObject::index_in_front_of_index(IndexTrio &one, IndexTrio &two)
{
	return (one.z < two.z);
}

vec3 vec_from_pos(GLfloat *pos)
{
	vec3 tmpVec = make_vec3(pos[0], pos[1],
	                        pos[2]);

	return tmpVec;
}


void GLObject::reorderIndices()
{
	_temp.resize(_indices.size() / 3);
	
	int count = 0;
	for (int i = 0; i < _indices.size(); i+=3)
	{
		int n = _indices[i];
		vec3 tmpVec = vec_from_pos(_vertices[n].pos);
		n = _indices[i + 1];
		vec3 tmpVec1 = vec_from_pos(_vertices[n].pos);
		n = _indices[i + 2];
		vec3 tmpVec2 = vec_from_pos(_vertices[n].pos);
		vec3_add_to_vec3(&tmpVec, tmpVec1);
		vec3_add_to_vec3(&tmpVec, tmpVec2);
		tmpVec = mat4x4_mult_vec(modelMat, tmpVec);
		_temp[count].z = tmpVec.z;
		_temp[count].index[0] = _indices[i];
		_temp[count].index[1] = _indices[i + 1];
		_temp[count].index[2] = _indices[i + 2];
		count++;
	}
	
	if (_backToFront)
	{
		std::sort(_temp.begin(), _temp.end(), index_behind_index);
	}
	else
	{
		std::sort(_temp.begin(), _temp.end(), index_in_front_of_index);
	}

	count = 0;

	for (int i = 0; i < _temp.size(); i++)
	{
		_indices[count + 0] = _temp[i].index[0];
		_indices[count + 1] = _temp[i].index[1];
		_indices[count + 2] = _temp[i].index[2];
		count += 3;
	}
}
