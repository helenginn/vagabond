#ifndef __vagabond_InkBond_fsh__
#define __vagabond_InkBond_fsh__

inline std::string InkBond_fsh()
{
	std::string str = 
	"varying vec4 vColor;\n"\
	"varying vec2 vTex;\n"\
	"varying vec4 vPos;\n"\
	"\n"\
	"uniform sampler2D bondTexture;\n"\
	"uniform vec3 focus;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	gl_FragColor = vColor;\n"\
	"	vec4 temp = texture2D(bondTexture, vTex);\n"\
	"	gl_FragColor = temp * vColor * 2.0;\n"\
	"	\n"\
	"	if (vColor[0] < 0.7 && vColor[1] < 0.7 && vColor[2] < 0.7) {\n"\
	"		gl_FragColor /= 2.0;\n"\
	"   }\n"\
	"	if (gl_FragColor[3] < 0.5) {\n"\
	"		discard;\n"\
	"	}\n"\
	"	float min_distance = -20.;\n"\
	"	float max_distance = -100.;\n"\
	"	if (focus[2] > -25.)\n"\
	"	{\n"\
	"		min_distance = focus[2] + 0.;\n"\
	"		max_distance = focus[2] - 8.;\n"\
	"	}\n"\
	"	if (vPos[2] > -2.) {\n"\
	"		discard;\n"\
	"	}\n"\
	"   float transparency = (vPos[2] - min_distance) / (max_distance - min_distance);\n"\
	"	transparency = max(transparency, 0.0);\n"\
	"	transparency = min(transparency, 1.0);\n"\
	"   gl_FragColor[3] = 1. - transparency;\n"\
	"	if (focus[2] < -120.)\n"\
	"	{\n"\
	"       gl_FragColor[3] = 1.;\n"\
	"   }\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return str;
}


#endif
