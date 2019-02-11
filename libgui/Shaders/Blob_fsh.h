#ifndef __vagabond_Blob_fsh__
#define __vagabond_Blob_fsh__

std::string Blob_fsh =
"#version 120\n"\
"varying vec4 vColor;\n"\
"varying vec2 vTex;\n"\
"\n"\
"uniform vec3 light_pos;\n"\
"\n"\
"uniform sampler2D blobTexture;\n"\
"\n"\
"void main()\n"\
"{\n"\
"	vec2 frag = gl_PointCoord;\n"\
"	vec2 xy = vec2(frag[0], frag[1]);\n"\
"	vec4 temp = texture2D(blobTexture, xy);\n"\
"	if (temp[3] < 0.05) {\n"\
"		discard;\n"\
"	}\n"\
"	gl_FragColor = temp * vColor;\n"\
"\n"\
"\n"\
"\n"\
"}\n";


#endif

