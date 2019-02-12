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
"uniform vec3 focus;\n"\
"\n"\
"void main()\n"\
"{\n"\
"   vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
"   gl_Position = projection * model * pos;\n"\
"   vec4 model4 = model * pos;\n"\
"	gl_PointSize = -500. / model4[2];\n"\
"	float min_distance = -20.;\n"\
"	float max_distance = -100.;\n"\
"	if (focus[2] > -15.)\n"\
"	{\n"\
"		min_distance = focus[2] + 0.;\n"\
"		max_distance = focus[2] - 8.;\n"\
"	}\n"\
"   float transparency = (model4[2] - min_distance) / (max_distance - min_distance);\n"\
"	transparency = max(transparency, 0.0);\n"\
"	transparency = min(transparency, 1.0);\n"\
"	vColor = color;\n"\
"   vColor[3] = 1. - transparency;\n"\
"}";


#endif

