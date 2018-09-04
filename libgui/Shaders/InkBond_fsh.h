#ifndef __vagabond_InkBond_fsh__
#define __vagabond_InkBond_fsh__

std::string InkBond_fsh =
"varying vec4 vColor;\n"\
"varying vec2 vTex;\n"\
"varying vec4 vPos;\n"\
"\n"\
"uniform sampler2D bondTexture;\n"\
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
"	if (vPos[2] < -15.) {\n"\
"		float dist = -vPos[2] - 15.;\n"\
"		float frac = 1. - min(dist / 6., 1.0);\n"\
"		gl_FragColor[3] *= frac;\n"\
"	}\n"\
"\n"\
"\n"\
"\n"\
"}\n";


#endif
