#ifndef __vagabond_HKL_vsh__
#define __vagabond_HKL_vsh__

inline std::string& hklVsh()
{
	static std::string HKL_vsh =
	"attribute vec3 normal;\n"\
	"attribute vec3 position;\n"\
	"attribute vec4 color;\n"\
	"\n"\
	"varying vec4 vColor;\n"\
	"varying vec4 vPos;\n"\
	"\n"\
	"uniform mat4 model;\n"\
	"uniform mat4 projection;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"   vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
	"   gl_Position = model * pos;\n"\
	"	gl_PointSize = 50. * color[3];\n"\
	"	vColor = color;\n"\
	"	vPos = model * pos;\n"\
	"}";
	return HKL_vsh;
}


#endif


