/**
 Frameworks for OpenGL

 */

#ifndef FRAMEWORKS_H

#include <QtGui/qopengl.h>
#include <QtGui/qopenglfunctions.h>
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
	GLfloat extra[4];
	GLfloat tex[2];
} Vertex;

typedef struct
{
	GLuint index[3];
	GLfloat z;
} IndexTrio;

class GLObject;
class Vagabond2GL;
class Density2GL;
typedef boost::shared_ptr<GLObject> GLObjectPtr;
typedef boost::shared_ptr<Vagabond2GL> Vagabond2GLPtr;
typedef boost::shared_ptr<Density2GL> Density2GLPtr;

#define ToDensity2GLPtr(a) (boost::static_pointer_cast<Density2GL>((a)))

#define FRAMEWORKS_H
#endif
