#ifndef __vagabond_Blob_vsh__
#define __vagabond_Blob_vsh__

std::string Blob_vsh =
"attribute vec3 normal;\n"\
"attribute vec3 position;\n"\
"attribute vec4 color;\n"\
"attribute vec2 tex;\n"\
"\n"\
"varying vec4 vColor;\n"\
"varying vec2 vTex;\n"\
"\n"\
"uniform mat4 projection;\n"\
"uniform mat4 model;\n"\
"uniform vec3 light_pos;\n"\
"\n"\
"void main()\n"\
"{\n"\
"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
"    gl_Position = projection * model * pos;\n"\
"    vec4 model4 = model * pos;\n"\
"	 gl_PointSize = -300. / model4[2];\n"\
"	 vColor = color;\n"\
"}";


#endif

