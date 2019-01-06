layout(location = 0) out vec3 delta_irradiance;
layout(location = 1) out vec3 irradiance;
uniform sampler3D single_rayleigh_scattering_texture;
uniform sampler3D single_mie_scattering_texture;
uniform sampler3D multiple_scattering_texture;
uniform int scattering_order;

void main() {
	delta_irradiance = ComputeIndirectIrradianceTexture(
      	ATMOSPHERE, single_rayleigh_scattering_texture,
      	single_mie_scattering_texture, multiple_scattering_texture,
      	gl_FragCoord.xy, scattering_order);
	irradiance = delta_irradiance;
}