#ifndef __vagabond_Blob_fsh__
#define __vagabond_Blob_fsh__

inline std::string& blobFsh()
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

