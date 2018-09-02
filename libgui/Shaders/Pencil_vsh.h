#ifndef __vagabond_Pencil_vsh__
#define __vagabond_Pencil_vsh__

std::string Pencil_vsh =
"attribute vec3 normal;\n"\
"attribute vec3 position;\n"\
"attribute vec4 color;\n"\
"attribute vec4 extra;\n"\
"attribute vec2 tex;\n"\
"\n"\
"varying vec4 vColor;\n"\
"varying vec2 vTex;\n"\
"\n"\
"uniform mat4 projection;\n"\
"uniform mat4 model;\n"\
"\n"\
"void main()\n"\
"{\n"\
"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
"    vec4 norm = vec4(normal[0], normal[1], normal[2], 1.0);\n"\
"    vec4 norm4 = model * norm;\n"\
"    float size = length(norm);\n"\
"    float redness = (5.0 - size) / 5.0;\n"\
"    if (redness < 0.) redness = 0.;\n"\
"    vec4 tmpColor = color;\n"\
"    tmpColor *= 1. - redness;"
"    tmpColor += redness * vec4(1., 0., 0., 1.);\n"\
"    vec4 modelPos = model * pos;\n"\
"    gl_Position = projection * modelPos;\n"\
"	 vColor = tmpColor;\n"\
"}";


#endif
