#ifndef __vagabond_Pencil_fsh__
#define __vagabond_Pencil_fsh__

std::string Pencil_fsh =
"varying vec4 vColor;\n"\
"varying vec2 vTex;\n"\
"varying vec3 vNormal;\n"\
"varying vec4 vPos;\n"\
"varying vec4 tPos;\n"\
"\n"\
"uniform sampler2D pencilTexture;\n"\
"uniform vec3 light_pos;\n"\
"\n"\
"void main()\n"\
"{\n"\
"	if (tPos[2] > -8.) {\n"\
"		discard;\n"\
"	}\n"\
"	if (tPos[2] < -20.) {\n"\
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
"	vec4 temp = texture2D(pencilTexture, tex);\n"\
"	gl_FragColor = temp;\n"\
"	gl_FragColor[0] *= vColor[0];\n"\
"	gl_FragColor[1] *= vColor[1];\n"\
"	gl_FragColor[2] *= vColor[2];\n"\
"	gl_FragColor[3] = 0.9 * diff;\n"\
"\n"\
"\n"\
"\n"\
"}\n";


#endif
