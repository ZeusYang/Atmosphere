#version 330
layout (location = 0) in vec3 position;

out vec3 view_ray;

uniform mat4 model_from_view;// inverses of viewMatrix
uniform mat4 view_from_clip;// inverses of projectMatrix

void main(){
	vec4 vertex = vec4(position,1.0f);
	view_ray = (model_from_view * vec4((view_from_clip * vertex).xyz, 0.0)).xyz;
	gl_Position = vertex;
}