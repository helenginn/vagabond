#ifndef __cluster_pointfsh__
#define __cluster_pointfsh__

inline std::string& pointFsh()
{
	static std::string Blob_fsh =
	"#version 330 core\n"\
	"in vec4 vColor;\n"\
	"in vec4 vPos;\n"\
	"\n"\
	"out vec4 fragColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	fragColor = vColor;\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return Blob_fsh;
}

#endif

