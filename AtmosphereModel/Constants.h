#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Atmosphere{

// 光学长度纹理size
constexpr int TRANSMITTANCE_TEXTURE_WIDTH   = 256;
constexpr int TRANSMITTANCE_TEXTURE_HEIGHT  = 64;

// 内散射(r,mu,mu_s,nu)的size
constexpr int SCATTERING_TEXTURE_R_SIZE     = 32;
constexpr int SCATTERING_TEXTURE_MU_SIZE    = 128;
constexpr int SCATTERING_TEXTURE_MU_S_SIZE  = 32;
constexpr int SCATTERING_TEXTURE_NU_SIZE    = 8;

// 内散射积分纹理size,3D纹理
constexpr int SCATTERING_TEXTURE_WIDTH      =
        SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;
constexpr int SCATTERING_TEXTURE_HEIGHT     = SCATTERING_TEXTURE_MU_SIZE;
constexpr int SCATTERING_TEXTURE_DEPTH      = SCATTERING_TEXTURE_R_SIZE;

// 辐照度纹理size
constexpr int IRRADIANCE_TEXTURE_WIDTH      = 64;
constexpr int IRRADIANCE_TEXTURE_HEIGHT     = 16;
}

#endif // CONSTANTS_H
