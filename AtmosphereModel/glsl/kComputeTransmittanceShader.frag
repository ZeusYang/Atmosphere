layout(location = 0) out vec3 transmittance;
void main() {/* 计算光学长度 */
	transmittance = ComputeTransmittanceToTopAtmosphereBoundaryTexture(
								ATMOSPHERE, gl_FragCoord.xy);
}