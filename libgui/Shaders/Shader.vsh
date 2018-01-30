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
    vec4 lightColor = vec4(1.0, 1.0, 1.0, 1.0);
    vec4 lightPos = vec4(0.0, 0.0, 0.0, 0.0);
    lightPos *= projection * model;
    vec4 lightVector = position - lightPos;
    vec4 unitLight = normalize(LightVector);
    float cosine = dot(unitLight, normal);
    float lambert = max(cosine, 0.0);
    float distance = distance(lightPos, position);
    float luminosity = 1 / (1 + distance * distance * 0.25);

    vec4 pos = vec4(position[0], position[1], position[2], 1.0);
    gl_Position = projection * model * pos;
    vColor = color * lightColor * lambert * luminosity;
}
