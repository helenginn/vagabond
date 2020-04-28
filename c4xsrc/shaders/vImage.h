#ifndef __Blot__Image_vsh__
#define __Blot__Image_vsh__

std::string vImage()
{
	std::string str =
	"attribute vec3 normal;\n"\
	"attribute vec3 position;\n"\
	"attribute vec4 color;\n"\
	"\n"\
	"uniform mat4 model;\n"\
	"\n"\
	"varying vec4 vColor;\n"\
	"varying vec4 vPos;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"    vec4 pos = model * vec4(position, 1.0);\n"\
	"    vPos = pos;\n"\
	"    gl_Position = pos;\n"\
	"	 vColor = color;\n"\
	"}";
	return str;
}


#endif
