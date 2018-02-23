//
//  shader.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 02/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

/*
 * Copyright (C) 2010 Josh A. Beam
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "shader.h"

#include <string.h>
#include <string>
#include "../libsrc/FileReader.h"

/*
 * Returns a string containing the text in
 * a vertex/fragment shader source file.
 */
static char *shaderLoadSource(const char *filePath)
{
    std::string contents = get_file_contents(filePath);
    return (char *)contents.c_str();
}

Shader::Shader()
{
    initializeOpenGLFunctions();
}

/*
 * Returns a shader object containing a shader
 * compiled from the given GLSL shader file.
 */
GLuint Shader::shaderCompileFromFile(GLenum type, const char *filePath, bool isString)
{
    GLuint shader;
    GLint length;
    GLint result;

    const char *cstr = NULL;
    
    if (!isString)
    {
        std::string contents = get_file_contents(filePath);
        cstr = contents.c_str();
        length = (GLint)contents.size();
    }
    else
    {
        cstr = filePath;
        length = (GLint)strlen(filePath);
    }

    /* create shader object, set the source, and compile */
    shader = glCreateShader(type);
    glShaderSource(shader, 1, &cstr, &length);
    glCompileShader(shader);

    /* make sure the compilation was successful */
    glGetShaderiv(shader, GL_COMPILE_STATUS, &result);
    if(result == GL_FALSE) {
        char *log;

        /* get the shader info log */
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
        log = (char *)malloc(length);
        glGetShaderInfoLog(shader, length, &result, log);

        /* print an error message and the info log */
        fprintf(stderr, "shaderCompileFromFile(): Unable to compile %s: %s\n", filePath, log);
        free(log);

        glDeleteShader(shader);
        return 0;
    }

    return shader;
}

/*
 * Compiles and attaches a shader of the
 * given type to the given program object.
 */
void Shader::shaderAttachFromFile(GLuint program, GLenum type, const char *filePath, bool isString)
{
    /* compile the shader */
    Shader me;
    GLuint shader = me.shaderCompileFromFile(type, filePath, isString);
    if(shader != 0) {
        /* attach the shader to the program */
        me.glAttachShader(program, shader);

        /* delete the shader - it won't actually be
         * destroyed until the program that it's attached
         * to has been destroyed */
        me.glDeleteShader(shader);
    }
}
