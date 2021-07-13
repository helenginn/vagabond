#ifndef __cluster__blobvsh__
#define __cluster__blobvsh__

inline std::string& pointVsh()
{
	static std::string Blob_vsh =
	"#version 330 core\n"\
	"in vec3 normal;\n"\
	"in vec3 position;\n"\
	"in vec4 color;\n"\
	"\n"\
	"out vec4 vColor;\n"\
	"out vec4 vPos;\n"\
	"out vec3 mPos;\n"\
	"\n"\
	"uniform mat4 model;\n"\
	"uniform mat4 projection;\n"\
	"uniform float size;\n"\
	"uniform int depth;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"   vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
	"   gl_Position = projection * model * pos;\n"\
	"	vPos = model * pos;\n"\
	"	mPos = mat3(model) * vec3(pos);\n"\
	"	float pointsize = size;\n"\
	"	if (depth > 0)\n"
	"	{\n"
	"		float frac = min(-mPos.z, 1.);\n"
	"		frac = max(frac, -1.);\n"
	"		frac = 0.5 - (frac / 3.5);\n"
	"		pointsize *= frac;\n"
	"	}\n"
	"	gl_PointSize = pointsize;\n"
	"	vColor = color;\n"\
	"}";
	return Blob_vsh;
}


#endif

