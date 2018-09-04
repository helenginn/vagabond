#ifndef __vagabond_InkBond__
#define __vagabond_InkBond__

std::string InkBond_vsh =
"attribute vec3 normal;\n"\
"attribute vec3 position;\n"\
"attribute vec4 color;\n"\
"attribute vec4 extra;\n"\
"attribute vec2 tex;\n"\
"\n"\
"varying vec4 vColor;\n"\
"varying vec2 vTex;\n"\
"varying vec4 vPos;\n"\
"\n"\
"uniform mat4 projection;\n"\
"uniform mat4 model;\n"\
"\n"\
"void main()\n"\
"{\n"\
"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
"    vec4 norm = vec4(normal[0], normal[1], normal[2], 1.0);\n"\

"	 vec4 norm4 = projection * model * norm;\n"\
"    vec4 model4 = projection * model * pos;\n"\
"    vec2 axis2 = vec2(norm4[0] - model4[0], norm4[1] - model4[1]);\n"\
"    vec2 axis = normalize(axis2);\n"\
"    vec2 clean_shift = vec2(-0.2, -0.3);\n"\
"    vColor = vec4(color[0], color[1], color[2], 1.0);\n"\
"    if (extra[0] < 0.5) {\n"\
"        clean_shift[1] *= -1.;\n"\
"    }\n"\
"    if (extra[1] > 0.5) {\n"\
"        clean_shift[0] = 0.0;\n"\
"    }\n"\
"    mat2 bondMat = mat2(axis[0], axis[1],\n"\
"	     	 			axis[1], -axis[0]);\n"\
"    vec2 shifted = bondMat * clean_shift;\n"\
"    model4[0] += shifted[0];\n"\
"    model4[1] += shifted[1];\n"\
"    gl_Position = model4;\n"\
"    vTex = tex;\n"\
"    vPos = model * pos;\n"\
"}";


#endif

 
