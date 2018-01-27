/**
 Frameworks for OpenGL

 */

#ifndef FRAMEWORKS_H

#include <OpenGL/OpenGL.h>
#include <OpenGL/gl.h>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#define MAX_LIGHTS 20
#define START_Z -64

typedef struct
{
	GLfloat pos[3];
	GLfloat normal[3];
	GLfloat color[4];
} Vertex;

class GLObject;
class Vagabond2GL;
typedef boost::shared_ptr<GLObject> GLObjectPtr;
typedef boost::shared_ptr<Vagabond2GL> Vagabond2GLPtr;

#define FRAMEWORKS_H
#endif
