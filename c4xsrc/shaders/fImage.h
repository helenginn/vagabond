#ifndef __Blot__Image_fsh__
#define __Blot__Image_fsh__

std::string fImage()
{
	std::string str =
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
	return str;
}


#endif
