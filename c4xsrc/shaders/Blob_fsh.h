#ifndef __cluster_pointfsh__
#define __cluster_pointfsh__

inline std::string& pointFsh()
{
	static std::string Blob_fsh =
	"#version 330 core\n"\
	"in vec4 vColor;\n"\
	"in vec4 vPos;\n"\
	"in vec3 mPos;\n"\
	"\n"\
	"uniform int depth;\n"\
	"\n"\
	"out vec4 fragColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"\n"\
	"	vec2 pos = gl_PointCoord - vec2(0.5, 0.5);\n"\
	"	if (length(pos) > 0.5) discard;\n"\
	"	if (vColor.a < 0.02) discard;\n"\
	"	fragColor = vColor;\n"\
	"\n"\
	"	if (false && depth > 0)\n"\
	"	{\n"\
	"		float amount = min((mPos.z + 1) / 2, 1);\n"\
	"		fragColor.xyz += (1 - fragColor.xyz) * amount;\n"\
	"\n"\
	"\n"\
	"	}\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return Blob_fsh;
}

#endif

