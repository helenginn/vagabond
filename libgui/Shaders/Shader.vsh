//
//  Shader.vsh
//  Vertex shader
//
//  Created by Helen Ginn on 18/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

attribute vec3 normal;
attribute vec3 position;
attribute vec4 color;

varying vec4 vColor;

uniform mat4 projection;
uniform mat4 model;

void main()
{
	vec4 pos = vec4(position[0], position[1], position[2], 1.0);
	gl_Position = projection * model * pos;
	vColor = color;
}
