inline std::string Shader_fsh()
{
	std::string str =
	"varying vec4 vColor;\n"\
	"varying vec4 vPos;\n"\
	"uniform vec3 focus;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	float min_distance = 200.;\n"\
	"	float max_distance = 400.;\n"\
	"	if (focus[2] < 100.)\n"\
	"	{\n"\
	"		min_distance = 70.;\n"\
	"		max_distance = 110.;\n"\
	"	}\n"\
	"	if (vPos[2] > max_distance) {\n"\
	"		discard;\n"\
	"	}\n"\
	"   float transparency = (vPos[2] - min_distance) / (max_distance - min_distance);\n"\
	"	transparency = max(transparency, 0.0);\n"\
	"	transparency = min(transparency, 1.0);\n"\
	"	vec4 tmpColor = vColor;\n"\
	"   for (int i = 0; i < 3; i++)\n"\
	"   {\n"\
	"       tmpColor[i] = tmpColor[i] + (1. - tmpColor[i]) * transparency;\n"\
	"   }\n"\
	"   tmpColor[3] = 0.85;\n"\
	"	gl_FragColor = tmpColor;\n"\
	"}\n";
	return str;
}
