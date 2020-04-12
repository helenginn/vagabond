#ifndef __Blot__Wipe_fsh__
#define __Blot__Wipe_fsh__

std::string fWipe =
"#version 120\n"\
"varying vec4 vColor;\n"\
"varying vec2 vTex;\n"\
"varying vec4 vPos;\n"\
"varying float vTime;\n"\
"\n"\
"uniform sampler2D pic_tex;\n"\
"\n"\
"void main()\n"\
"{\n"\
"	vec4 frag = gl_FragCoord;\n"\
"	vec2 tex = vec2(vTex[0], vTex[1]);\n"\
"	vec4 temp = texture2D(pic_tex, tex);\n"\
"	gl_FragColor = vec4(temp[2], temp[1], temp[0], temp[3]);\n"\
"	gl_FragColor += vColor;\n"\
"	if (frag[1] > vTime)\n"\
"	{\n"\
"		gl_FragColor[3] = 0.0;\n"\
"	}\n"\
"\n"\
"\n"\
"\n"\
"}\n";


#endif
