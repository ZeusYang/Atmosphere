uniform sampler2D transmittance_texture;
uniform sampler3D scattering_texture;
uniform sampler3D single_mie_scattering_texture;
uniform sampler2D irradiance_texture;

RadianceSpectrum GetSolarRadiance() {
	return ATMOSPHERE.solar_irradiance /
      		(PI * ATMOSPHERE.sun_angular_radius * ATMOSPHERE.sun_angular_radius);
}
RadianceSpectrum GetSkyRadiance(
	Position camera, Direction view_ray, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance) {
	return GetSkyRadiance(ATMOSPHERE, transmittance_texture,
      		scattering_texture, single_mie_scattering_texture,
      		camera, view_ray, shadow_length, sun_direction, transmittance);
}
RadianceSpectrum GetSkyRadianceToPoint(
	Position camera, Position point, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance) {
	return GetSkyRadianceToPoint(ATMOSPHERE, transmittance_texture,
      		scattering_texture, single_mie_scattering_texture,
      		camera, point, shadow_length, sun_direction, transmittance);
}
IrradianceSpectrum GetSunAndSkyIrradiance(
	Position p, Direction normal, Direction sun_direction,
	out IrradianceSpectrum sky_irradiance) {
	return GetSunAndSkyIrradiance(ATMOSPHERE, transmittance_texture,
      		irradiance_texture, p, normal, sun_direction, sky_irradiance);
}
