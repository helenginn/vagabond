#ifndef __cluster_pointfsh__
#define __cluster_pointfsh__

inline std::string& pointFsh()
{
	static std::string Blob_fsh =
	"#version 120\n"\
	"varying vec4 vColor;\n"\
	"varying vec4 vPos;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	gl_FragColor = vColor;\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return Blob_fsh;
}

#endif

