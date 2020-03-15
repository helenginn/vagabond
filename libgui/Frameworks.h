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
#include "../libsrc/vec3.h"

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

inline void vec3ToVertex(Vertex &v, vec3 &vec)
{
	v.pos[0] = vec.x;
	v.pos[1] = vec.y;
	v.pos[2] = vec.z;
}

inline void vec3ToNormal(Vertex &v, vec3 &vec)
{
	v.normal[0] = vec.x;
	v.normal[1] = vec.y;
	v.normal[2] = vec.z;
}

typedef struct
{
	GLuint index[3];
	GLfloat z;
} IndexTrio;

class GLObject;
class Vagabond2GL;
class Density2GL;
class Bonds2GL;
class Atoms2GL;
class Multi2GL;
class Connect2GL;
class Selected2GL;
typedef boost::shared_ptr<GLObject> GLObjectPtr;
typedef boost::shared_ptr<Vagabond2GL> Vagabond2GLPtr;
typedef boost::shared_ptr<Bonds2GL> Bonds2GLPtr;
typedef boost::shared_ptr<Atoms2GL> Atoms2GLPtr;
typedef boost::shared_ptr<Multi2GL> Multi2GLPtr;
typedef boost::shared_ptr<Connect2GL> Connect2GLPtr;
typedef boost::shared_ptr<Selected2GL> Selected2GLPtr;
typedef boost::shared_ptr<Density2GL> Density2GLPtr;

#define ToVagabond2GLPtr(a) (boost::static_pointer_cast<Vagabond2GL>((a)))
#define ToDensity2GLPtr(a) (boost::static_pointer_cast<Density2GL>((a)))

#define FRAMEWORKS_H
#endif
