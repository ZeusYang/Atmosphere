#version 330
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 texcoord;


uniform mat4 modelMatrix;
uniform mat4 viewProjectMatrix;

out vec3 Normal;

void main(){
	Normal = normal;
	gl_Position = viewProjectMatrix*modelMatrix*vec4(position, 1.0f);
}