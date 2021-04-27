#ifndef __vagabond_Pencil_fsh__
#define __vagabond_Pencil_fsh__

inline std::string Pencil_fsh()
{
	std::string str = 
"\n"\
"#version 330 core\n"\
"in vec4 vColor;\n"\
"in vec2 vTex;\n"\
"in vec3 vNormal;\n"\
"in vec4 vExtra;\n"\
"in vec4 vPos;\n"\
"in vec4 tPos;\n"\
"\n"\
"uniform sampler2D pencilTexture;\n"\
"uniform vec3 light_pos;\n"\
"uniform vec3 focus;\n"\
"\n"\
"out vec4 fragColor;\n"\
"\n"\
"\n"\
"void main()\n"\
"{\n"\
"	if (tPos[2] > -4.) {\n"\
"		discard;\n"\
"	}\n"\
"	vec3 pos = vec3(vPos[0], vPos[1], vPos[2]);\n"\
"	vec3 lightDir = normalize(pos - light_pos);\n"\
"	float diff = abs(dot(vNormal, lightDir));\n"\
"	diff = sqrt(1. - diff * diff);\n"\
"   float cutoff = 0.4;\n"\
"   if (diff < cutoff) {\n"\
"       diff = 0.;\n"\
"	} else {\n"\
"	    diff -= cutoff; diff /= 1. - cutoff;\n"\
"	}\n"\
"	vec2 tex = vec2(tPos[0], tPos[1]);\n"\
"	vec4 temp = texture(pencilTexture, tex);\n"\
"	fragColor = temp;\n"\
"	fragColor[0] *= vColor[0];\n"\
"	fragColor[1] *= vColor[1];\n"\
"	fragColor[2] *= vColor[2];\n"\
"	fragColor[3] = 0.9 * diff;\n"\
"	float min_distance = -focus[2] + 5.;\n"\
"	float max_distance = -focus[2] - 5.;\n"\
"	if (focus[2] < +15.)\n"\
"	{\n"\
"		min_distance = -focus[2] + 2.;\n"\
"		max_distance = -focus[2] - 3.;\n"\
"	}\n"\
"	if ((tPos[2] < max_distance || tPos[2] > min_distance)) \n"\
"	{ discard; } \n"\
"\n"\
"\n"\
"\n"\
"\n"\
"}\n";
return str;
}


#endif
