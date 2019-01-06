#version 330

in vec3 Normal;
out vec4 color;
void main(){
	color = vec4(1.0,0.0,0.0,1.0);
	color = vec4(Normal,1.0f);
}