#ifndef ATMOSPHEREMODEL_H
#define ATMOSPHEREMODEL_H

#include <string>
#include <vector>
#include <array>
#include <functional>
#include <QOpenGLShader>
#include <QOpenGLFunctions>
#include "Utilities/Entity.h"
#include "Utilities/OpenGLTexture.h"
#include "Utilities/ShaderProgram.h"

namespace Atmosphere{

/**********************************
 * 大气层，宽度为width，密度定义为
 * 'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term'
 * 【0，1】，其中h为海拔高度
 **********************************/
class DensityProfileLayer {
public:
    DensityProfileLayer() : DensityProfileLayer(0.0, 0.0, 0.0, 0.0, 0.0) {}
    DensityProfileLayer(double width, double exp_term, double exp_scale,
                        double linear_term, double constant_term)
        : width(width), exp_term(exp_term), exp_scale(exp_scale),
          linear_term(linear_term), constant_term(constant_term) {}

    // 海拔高度h处的空气分子(气溶胶)密度=
    // 'exp_term' * exp('exp_scale' * h)
    //+ 'linear_term' * h
    //+ 'constant_term'
    double width;//大气层宽度
    double exp_term;
    double exp_scale;
    double linear_term;
    double constant_term;
};

class AtmosphereModel {
public:
    AtmosphereModel(
            // 太阳波长，单位nm
            const std::vector<double>& wavelengths,
            // 太阳辐照度
            const std::vector<double>& solar_irradiance,
            // 太阳角半径
            double sun_angular_radius,
            // 大气层底层到星球中心的距离(内半径)
            double bottom_radius,
            // 大气层外层到星球中心的距离(外半径)
            double top_radius,
            // 大气空气分子密度分布
            const std::vector<DensityProfileLayer>& rayleigh_density,
            // (空气分子)rayleigh散射系数=海拔高度h处的rayleigh_scattering*rayleigh_density
            const std::vector<double>& rayleigh_scattering,
            // 大气气溶胶密度分布
            const std::vector<DensityProfileLayer>& mie_density,
            // (气溶胶)mie散射系数=海拔高度h处的mie_scattering*mie_density
            const std::vector<double>& mie_scattering,
            // (气溶胶)mie消光系数=海拔高度h处mie_extinction*mie_density
            const std::vector<double>& mie_extinction,
            // 气溶胶Cornette-Shanks的相位函数参数值g
            double mie_phase_function_g,
            // 大气中吸收光线的空气分子密度
            const std::vector<DensityProfileLayer>& absorption_density,
            // 大气中吸收光线的消光系数=海拔高度h处absorption_extinction*absorption_density
            const std::vector<double>& absorption_extinction,
            // 地面的平均反照率
            const std::vector<double>& ground_albedo,
            // 太阳最大的天顶角,弧度制
            double max_sun_zenith_angle,
            // 长度单位
            double length_unit_in_meters,
            // 是否把单次mie散射和rayleigh散射以及多次散射合并到一个纹理
            bool combine_scattering_textures,
            // 采用半精度浮点数还是单精度浮点数
            bool half_precision);

    ~AtmosphereModel();

    QOpenGLShader* getShader(){return atmosphere_shader_;}

    void init(unsigned int num_scattering_orders = 8);

    void setProgramUniforms(
            GLuint program,                                         //着色器程序
            GLuint transmittance_texture_unit,                             //透射率纹理单元
            GLuint scattering_texture_unit,                                //散射纹理单元
            GLuint irradiance_texture_unit,                                //辐照度纹理单元
            GLuint optional_single_mie_scattering_texture_unit = 0) const; //可选的单次mie散射纹理单元

    static constexpr double kLambdaR = 680.0;
    static constexpr double kLambdaG = 550.0;
    static constexpr double kLambdaB = 440.0;

private:
    typedef std::array<double, 3> vec3;
    typedef std::array<float, 9> mat3;

    void precompute(
            OpenGLTexture &delta_irradiance_texture,
            OpenGLTexture &delta_rayleigh_scattering_texture,
            OpenGLTexture &delta_mie_scattering_texture,
            OpenGLTexture &delta_scattering_density_texture,
            OpenGLTexture &delta_multiple_scattering_texture,
            const vec3& lambdas,
            unsigned int num_scattering_orders);

    bool half_precision_;
    bool rgb_format_supported_;
    std::function<std::string(const vec3&)> glsl_header_factory_;
    OpenGLTexture *transmittance_texture_;
    OpenGLTexture *scattering_texture_;
    OpenGLTexture *optional_single_mie_scattering_texture_;
    OpenGLTexture *irradiance_texture_;
    QOpenGLShader *atmosphere_shader_;

    Entity *quad;

    void drawQuad(const std::vector<bool>& enable_blend);

    double Interpolate(const std::vector<double>& wavelengths,
                        const std::vector<double>& wavelength_function,
                        double wavelength);
    bool IsFramebufferRgbFormatSupported(bool half_precision);

    QString getShaderFromFile(QString path);
};

}

#endif // ATMOSPHEREMODEL_H
