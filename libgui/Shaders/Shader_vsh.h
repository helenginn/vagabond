inline std::string Shader_vsh()
{
	std::string str =
	"attribute vec3 normal;\n"\
	"attribute vec3 position;\n"\
	"attribute vec4 color;\n"\
	"\n"\
	"varying vec4 vColor;\n"\
	"varying vec4 vPos;\n"\
	"\n"\
	"uniform mat4 projection;\n"\
	"uniform mat4 model;\n"\
	"uniform vec3 focus;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"   mat3 normMat = mat3(model[0][0], model[0][1], model[0][2],\n"\
	"                       model[1][0], model[1][1], model[2][2],\n"\
	"                       model[2][0], model[2][1], model[2][2]);\n"\
	"   vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
	"   vec4 modelPos = model * pos;\n"\
	"   vec3 norm3 = normalize(normMat * normal);\n"\
	"   vec4 norm = vec4(norm3[0], norm3[1], norm3[2], 1.0);\n"\
	"   gl_Position = projection * modelPos;\n"\
	"   vec4 lightColor = vec4(1.0, 1.0, 1.0, 1.0);\n"\
	"   vec4 lightPos = vec4(0.0, 0.0, 0.0, 0.0);\n"\
	"   vec4 lightVector = modelPos - lightPos;\n"\
	"   vec4 unitLight = normalize(lightVector);\n"\
	"   float cosine = 1.0 - abs(dot(unitLight, norm));\n"\
	"   float lambert = max(cosine, 0.0);\n"\
	"   float distance = gl_Position[2];\n"\
	"   float luminosity = 1.0 / (1. + distance * distance * 0.25);\n"\
	"   float complete = lambert * luminosity;\n"\
	"   complete *= complete;\n"\
	"   complete *= 0.3;\n"\
	"	vColor = color + complete * lightColor;\n"\
	"	vPos = modelPos;\n"\
	"}";
	return str;
}
