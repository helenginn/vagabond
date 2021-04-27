#ifndef __vagabond_Pencil_vsh__
#define __vagabond_Pencil_vsh__

inline std::string Pencil_vsh()
{
	std::string str = 
	"#version 330 core\n"\
	"in vec3 normal;\n"\
	"in vec3 position;\n"\
	"in vec4 color;\n"\
	"in vec4 extra;\n"\
	"in vec2 tex;\n"\
	"\n"\
	"out vec4 vColor;\n"\
	"out vec3 vNormal;\n"\
	"out vec4 vPos;\n"\
	"out vec4 tPos;\n"\
	"out vec2 vTex;\n"\
	"out vec4 vExtra;\n"\
	"\n"\
	"uniform mat4 projection;\n"\
	"uniform mat4 model;\n"\
	"uniform vec3 light_pos;\n"\
	"uniform vec3 focus;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
	"    vPos = pos;\n"\
	"    tPos = model * pos;\n"\
	"    gl_Position = projection * tPos;\n"\
	"    vNormal = normalize(normal);\n"\
	"	 vColor = color;\n"\
	"	 vExtra = extra;\n"\
	"	 vTex = vec2(tex[0], tex[1]);\n"\
	"}";
	return str;
}


#endif
