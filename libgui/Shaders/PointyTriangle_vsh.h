#ifndef __vagabond_PointyTriangle_vsh__
#define __vagabond_PointyTriangle_vsh__

inline std::string PointyTriangle_vsh ()
{
	std::string str =
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
"    float l = length(normal);\n"\
"    mat3 rot_only = mat3(1);\n"\
"    rot_only[0] = vec3(model[0][0], model[0][1], model[0][2]);\n"\
"    rot_only[1] = vec3(model[1][0], model[1][1], model[1][2]);\n"\
"    rot_only[2] = vec3(model[2][0], model[2][1], model[2][2]);\n"\
"    vec3 norm4 = rot_only * normal;\n"\
"    vec2 axis = normalize(vec2(norm4[0], norm4[1]));\n"\
"    vec2 clean_shift = vec2(extra[0], extra[1]);\n"\
"    mat2 dirMat = mat2(axis[0], axis[1],\n"\
"	     	 			-axis[1], axis[0]);\n"\
"    vec2 shifted = l * dirMat * clean_shift;\n"\
"    vec4 model4 = model * pos;\n"\
"    model4[0] += shifted[0];\n"\
"    model4[1] += shifted[1];\n"\
"    model4 = projection * model4;\n"\
"    gl_Position = model4;\n"\
"    vColor = color;\n"\
"    vTex = tex;\n"\
"    vPos = model * pos;\n"\
"}";
return str;
}


#endif

 


