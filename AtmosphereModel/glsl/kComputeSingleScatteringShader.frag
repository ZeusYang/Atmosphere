layout(location = 0) out vec3 delta_rayleigh;
layout(location = 1) out vec3 delta_mie;
layout(location = 2) out vec4 scattering;
layout(location = 3) out vec3 single_mie_scattering;
uniform sampler2D transmittance_texture;
uniform int layer;

void main() {
	ComputeSingleScatteringTexture(
      		ATMOSPHERE, transmittance_texture, vec3(gl_FragCoord.xy, layer+0.5f),
      		delta_rayleigh, delta_mie);
	scattering = vec4(delta_rayleigh.rgb,delta_mie.r);// rayleigh散射,alpha设为delta_mie.r是合并纹理的情况
	single_mie_scattering = delta_mie;// mie散射
}