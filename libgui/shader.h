//
//  shader.hpp
//  VagabondViewer
//
//  Created by Helen Ginn on 02/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef shader_hpp
#define shader_hpp

#include "Frameworks.h"


class Shader : public QObject, public QOpenGLFunctions
{
public:
    Shader();
    static void shaderAttachFromFile(GLuint program, GLenum type, const char *filePath, bool isString);

private:
	
    GLuint shaderCompileFromFile(GLenum type, const char *filePath, bool isString);
};

#endif /* shader_hpp */
