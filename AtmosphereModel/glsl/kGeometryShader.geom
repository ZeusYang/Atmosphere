#version 330
layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;
uniform int layer;
// 对于3D纹理，我们要确定写入的是哪一层
void main() {
	gl_Position = gl_in[0].gl_Position;
	gl_Layer = layer;
	EmitVertex();
	gl_Position = gl_in[1].gl_Position;
	gl_Layer = layer;
	EmitVertex();
	gl_Position = gl_in[2].gl_Position;
	gl_Layer = layer;
	EmitVertex();
	EndPrimitive();
}