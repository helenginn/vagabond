std::string Shader_vsh =
"attribute vec3 normal;\n"\
"attribute vec3 position;\n"\
"attribute vec4 color;\n"\
"\n"\
"varying vec4 vColor;\n"\
"\n"\
"uniform mat4 projection;\n"\
"uniform mat4 model;\n"\
"\n"\
"void main()\n"\
"{\n"\
"    mat3 normMat = mat3(model[0][0], model[0][1], model[0][2],\n"\
"                        model[1][0], model[1][1], model[2][2],\n"\
"                        model[2][0], model[2][1], model[2][2]);\n"\
"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
"    vec3 norm3 = normalize(normMat * normal);\n"\
"    vec4 norm = vec4(norm3[0], norm3[1], norm3[2], 1.0);\n"\
"    vec4 modelPos = model * pos;\n"\
"    gl_Position = projection * modelPos;\n"\
"    vec4 lightColor = vec4(1.0, 1.0, 1.0, 1.0);\n"\
"    vec4 lightPos = vec4(0.0, 0.0, 0.0, 0.0);\n"\
"    vec4 lightVector = modelPos - lightPos;\n"\
"    vec4 unitLight = normalize(lightVector);\n"\
"    float cosine = 1.0 - abs(dot(unitLight, norm));\n"\
"    float lambert = max(cosine, 0.0);\n"\
"    float distance = gl_Position[2];\n"\
"    float luminosity = 1.0 / (1. + distance * distance * 0.25);\n"\
"    float complete = lambert * luminosity;\n"\
"    complete *= complete;\n"\
"    complete *= 0.3;\n"\
"	 float min_distance = -20.;\n"\
"	 float max_distance = -100.;\n"\
"    float transparency = (modelPos[2] - min_distance) / (max_distance - min_distance);\n"\
"	 transparency = max(transparency, 0.0);\n"\
"	 transparency = min(transparency, 1.0);\n"\
"    vColor = color + complete * lightColor;\n"\
""\
"    for (int i = 0; i < 3; i++)\n"\
"    {\n"\
"        vColor[i] = vColor[i] + (1. - vColor[i]) * transparency;\n"\
"    }\n"\
"	 vColor[3] = 0.85;\n"\
"}";



 
