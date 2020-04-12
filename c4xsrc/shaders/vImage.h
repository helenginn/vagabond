#ifndef __Blot__Image_vsh__
#define __Blot__Image_vsh__

std::string vImage =
"attribute vec3 normal;\n"\
"attribute vec3 position;\n"\
"attribute vec4 color;\n"\
"attribute vec4 extra;\n"\
"attribute vec2 tex;\n"\
"\n"\
"uniform mat4 model;\n"\
"\n"\
"varying vec4 vColor;\n"\
"varying vec4 vPos;\n"\
"varying vec2 vTex;\n"\
"\n"\
"void main()\n"\
"{\n"\
"    vec4 pos = model * vec4(position, 1.0);\n"\
"    vPos = pos;\n"\
"    gl_Position = pos;\n"\
"	 vColor = color;\n"\
"	 vTex = tex;\n"\
"}";


#endif
