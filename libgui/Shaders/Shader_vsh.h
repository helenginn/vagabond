std::string Shader_vsh =
"attribute vec3 normal;\n"\
"attribute vec3 position;\n"\
"attribute vec4 color;\n"\
"\n"\
"varying vec4 vColor;\n"\
"\n"\
"uniform mat4 projection;\n"\
"uniform mat4 model;\n"\
"\n"\
"void main()\n"\
"{\n"\
"	vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
"	gl_Position = projection * model * pos;\n"\
"	vColor = color;\n"\
"}";

