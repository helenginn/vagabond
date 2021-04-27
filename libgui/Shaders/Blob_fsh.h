#ifndef __vagabond_Blob_fsh__
#define __vagabond_Blob_fsh__

inline std::string Blob_fsh()
{
	std::string str =
	"#version 330 core\n"\
	"in vec4 vColor;\n"\
	"in vec2 vTex;\n"\
	"in vec4 vPos;\n"\
	"\n"\
	"uniform vec3 light_pos;\n"\
	"uniform vec3 focus;\n"\
	"\n"\
	"uniform sampler2D blobTexture;\n"\
	"\n"\
	"out vec4 fragColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	vec2 frag = gl_PointCoord;\n"\
	"	vec2 xy = vec2(frag[0], frag[1]);\n"\
	"	vec4 temp = texture(blobTexture, xy);\n"\
	"	if (temp[3] < 0.05) {\n"\
	"		discard;\n"\
	"	}\n"\
	"	float min_distance = 200.;\n"\
	"	float max_distance = 400.;\n"\
	"	if (focus[2] < 100.)\n"\
	"	{\n"\
	"		min_distance = focus[2];\n"\
	"		max_distance = focus[2] + 5.;\n"\
	"	}\n"\
	"	if (vPos[2] > max_distance) {\n"\
	"		discard;\n"\
	"	}\n"\
	"   float transparency = (vPos[2] - min_distance) / (max_distance - min_distance);\n"\
	"	transparency = max(transparency, 0.3);\n"\
	"	transparency = min(transparency, 1.0);\n"\
	"	vec4 tmpColor = vColor;\n"\
	"   for (int i = 0; i < 3; i++)\n"\
	"   {\n"\
	"       tmpColor[i] = tmpColor[i] + (1. - tmpColor[i]) * transparency;\n"\
	"   }\n"\
	"   tmpColor[3] = 0.85;\n"\
	"	fragColor = tmpColor;\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";

	return str;
}


#endif

