// 以下define便于理解数字背后的含义
#define Length float
#define Wavelength	float
#define Angle	 float
#define SolidAngle float
#define Power	 float
#define LuminousPower float

#define Number float
#define InverseLength float
#define Area float
#define Volume float
#define NumberDensity float
#define Irradiance float
#define Radiance float
#define SpectralPower float
#define SpectralIrradiance float
#define SpectralRadiance float
#define SpectralRadianceDensity float
#define ScatteringCoefficient float
#define InverseSolidAngle float

// A generic function from Wavelength to some other type.
#define AbstractSpectrum vec3
// A function from Wavelength to Number.
#define DimensionlessSpectrum vec3
// A function from Wavelength to SpectralPower.
#define PowerSpectrum vec3
// A function from Wavelength to SpectralIrradiance.
#define IrradianceSpectrum vec3
// A function from Wavelength to SpectralRadiance.
#define RadianceSpectrum vec3
// A function from Wavelength to SpectralRadianceDensity.
#define RadianceDensitySpectrum vec3
// A function from Wavelength to ScaterringCoefficient.
#define ScatteringSpectrum vec3

// A position in 3D (3 length values).
#define Position vec3
// A unit direction vector in 3D (3 unitless values).
#define Direction vec3

#define TransmittanceTexture sampler2D
#define AbstractScatteringTexture sampler3D
#define ReducedScatteringTexture sampler3D
#define ScatteringTexture sampler3D
#define ScatteringDensityTexture sampler3D
#define IrradianceTexture sampler2D

const Length m = 1.0;//米
const Wavelength nm = 1.0;//纳米
const Angle rad = 1.0;//弧度
const SolidAngle sr = 1.0;//立体弧度
const Power watt = 1.0;//瓦特

const float PI = 3.14159265358979323846;

const Length km = 1000.0 * m;
const Area m2 = m * m;
const Volume m3 = m * m * m;
const Angle pi = PI * rad;
const Angle deg = pi / 180.0;
const Irradiance watt_per_square_meter = watt / m2;
const Radiance watt_per_square_meter_per_sr = watt / (m2 * sr);
const SpectralIrradiance watt_per_square_meter_per_nm = watt / (m2 * nm);
const SpectralRadiance watt_per_square_meter_per_sr_per_nm =
    watt / (m2 * sr * nm);
const SpectralRadianceDensity watt_per_cubic_meter_per_sr_per_nm =
    watt / (m3 * sr * nm);

struct DensityProfileLayer {//大气层密度剖面
	Length width;
	Number exp_term;
	InverseLength exp_scale;
	InverseLength linear_term;
	Number constant_term;
};

struct DensityProfile {
	DensityProfileLayer layers[2];
};

struct AtmosphereParameters {//大气层参数模型
	// 大气层顶部的太阳辐照度
	IrradianceSpectrum solar_irradiance;
	// 太阳角半径
  	Angle sun_angular_radius;
  	// 星球中心到大气层底部的距离,即地球半径
  	Length bottom_radius;
  	// 星球中心到大气层顶部的距离
  	Length top_radius;
  	// 空气分子的密度分布,[0,1]
  	DensityProfile rayleigh_density;
  	// 海拔为h处的rayleigh散射系数 = rayleigh_scattering * rayleigh_density
  	ScatteringSpectrum rayleigh_scattering;
  	// 气溶胶的密度剖面,[0,1]
  	DensityProfile mie_density;
	// mie散射系数 = mie_scattering * mie_density
  	ScatteringSpectrum mie_scattering;
  	// 气溶胶的消光系数 = mie_extinction * mie_density
  	ScatteringSpectrum mie_extinction;
  	// 气溶胶的Cornette-Shanks相位函数中的非对称参数
  	Number mie_phase_function_g;
  	// 吸收光线的空气分子的密度分布,[0,1]
  	DensityProfile absorption_density;
  	// 吸收光线的空气消光系数 = absorption_extinction * absorption_density
  	ScatteringSpectrum absorption_extinction;
  	// 地面的平均反照率(反射系数)
  	DimensionlessSpectrum ground_albedo;
  	// 太阳天最大天顶角的cos值(cos因此最小)，用于后面大气散射的预计算
  	Number mu_s_min;
};