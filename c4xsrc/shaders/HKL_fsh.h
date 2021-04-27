#ifndef __vagabond_HKL_fsh__
#define __vagabond_HKL_fsh__

inline std::string& hklFsh()
{
	static std::string HKL_fsh =
	"#version 330 core\n"\
	"in vec4 vColor;\n"\
	"in vec4 vPos;\n"\
	"\n"\
	"out vec4 fragColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	if (abs(vPos[2] / vPos[3]) > 0.05) {\n"\
	"		discard;\n"\
	"	}\n"\
	"	fragColor = vColor;\n"\
	"	vec2 frag = gl_PointCoord;\n"\
	"	vec2 xy = vec2(frag[0] - 0.5, frag[1] - 0.5);\n"\
	"	float dist = length(xy);\n"\
	"	if (dist > 0.5) {"
	"		discard;\n"\
	"	}\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return HKL_fsh;
}

#endif


