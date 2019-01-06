#include "AtmosphereModel.h"
#include "Constants.h"
#include <QFile>
#include <assert.h>
#include <QDebug>
#include <iostream>
#include "Utilities/Geometry/Quad.h"
#include "Utilities/OpenGLVersion.h"

#define F(c) c *f = (GlobalContext::contextFunc)

Atmosphere::AtmosphereModel::AtmosphereModel(
        const std::vector<double> &wavelengths,
        const std::vector<double> &solar_irradiance,
        double sun_angular_radius,
        double bottom_radius,
        double top_radius,
        const std::vector<Atmosphere::DensityProfileLayer> &rayleigh_density,
        const std::vector<double> &rayleigh_scattering,
        const std::vector<Atmosphere::DensityProfileLayer> &mie_density,
        const std::vector<double> &mie_scattering,
        const std::vector<double> &mie_extinction,
        double mie_phase_function_g,
        const std::vector<Atmosphere::DensityProfileLayer> &absorption_density,
        const std::vector<double> &absorption_extinction,
        const std::vector<double> &ground_albedo,
        double max_sun_zenith_angle,
        double length_unit_in_meters,
        bool combine_scattering_textures,
        bool half_precision):
    half_precision_(half_precision),
    rgb_format_supported_(IsFramebufferRgbFormatSupported(half_precision)),
    atmosphere_shader_(nullptr)
{
    auto to_string = [&](const std::vector<double>& v,
            const vec3& lambdas, double scale) {
        double r = Interpolate(wavelengths, v, lambdas[0]) * scale;
        double g = Interpolate(wavelengths, v, lambdas[1]) * scale;
        double b = Interpolate(wavelengths, v, lambdas[2]) * scale;
        return "vec3(" + std::to_string(r) + "," + std::to_string(g) + "," +
                std::to_string(b) + ")";
    };

    auto density_layer =
            [length_unit_in_meters](const DensityProfileLayer& layer) {
        return "DensityProfileLayer(" +
                std::to_string(layer.width / length_unit_in_meters) + "," +
                std::to_string(layer.exp_term) + "," +
                std::to_string(layer.exp_scale * length_unit_in_meters) + "," +
                std::to_string(layer.linear_term * length_unit_in_meters) + "," +
                std::to_string(layer.constant_term) + ")";
    };

    auto density_profile =
            [density_layer](std::vector<DensityProfileLayer> layers) {
        constexpr int kLayerCount = 2;
        while (layers.size() < kLayerCount) {
            layers.insert(layers.begin(), DensityProfileLayer());
        }
        std::string result = "DensityProfile(DensityProfileLayer[" +
                std::to_string(kLayerCount) + "](";
        for (int i = 0; i < kLayerCount; ++i) {
            result += density_layer(layers[i]);
            result += i < kLayerCount - 1 ? "," : "))";
        }
        return result;
    };
    // Shader文件的首部
    std::string definitions_glsl =
            getShaderFromFile(QString(":/AtmosphereModel/glsl/definitions.frag"))
            .toStdString();
    std::string functions_glsl =
            getShaderFromFile(QString(":/AtmosphereModel/glsl/functions.frag"))
            .toStdString();

    glsl_header_factory_ = [=](const vec3& lambdas) {
        return
                "#version 330\n"
                "#define IN(x) const in x\n"
                "#define OUT(x) out x\n"
                "const int TRANSMITTANCE_TEXTURE_WIDTH = " +
                std::to_string(TRANSMITTANCE_TEXTURE_WIDTH) + ";\n" +
                "const int TRANSMITTANCE_TEXTURE_HEIGHT = " +
                std::to_string(TRANSMITTANCE_TEXTURE_HEIGHT) + ";\n" +
                "const int SCATTERING_TEXTURE_R_SIZE = " +
                std::to_string(SCATTERING_TEXTURE_R_SIZE) + ";\n" +
                "const int SCATTERING_TEXTURE_MU_SIZE = " +
                std::to_string(SCATTERING_TEXTURE_MU_SIZE) + ";\n" +
                "const int SCATTERING_TEXTURE_MU_S_SIZE = " +
                std::to_string(SCATTERING_TEXTURE_MU_S_SIZE) + ";\n" +
                "const int SCATTERING_TEXTURE_NU_SIZE = " +
                std::to_string(SCATTERING_TEXTURE_NU_SIZE) + ";\n" +
                "const int IRRADIANCE_TEXTURE_WIDTH = " +
                std::to_string(IRRADIANCE_TEXTURE_WIDTH) + ";\n" +
                "const int IRRADIANCE_TEXTURE_HEIGHT = " +
                std::to_string(IRRADIANCE_TEXTURE_HEIGHT) + ";\n" +
                (combine_scattering_textures ?
                     "#define COMBINED_SCATTERING_TEXTURES\n" : "") +
                definitions_glsl +
                "const AtmosphereParameters ATMOSPHERE = AtmosphereParameters(\n" +
                to_string(solar_irradiance, lambdas, 1.0) + ",\n" +
                std::to_string(sun_angular_radius) + ",\n" +
                std::to_string(bottom_radius / length_unit_in_meters) + ",\n" +
                std::to_string(top_radius / length_unit_in_meters) + ",\n" +
                density_profile(rayleigh_density) + ",\n" +
                to_string(
                    rayleigh_scattering, lambdas, length_unit_in_meters) + ",\n" +
                density_profile(mie_density) + ",\n" +
                to_string(mie_scattering, lambdas, length_unit_in_meters) + ",\n" +
                to_string(mie_extinction, lambdas, length_unit_in_meters) + ",\n" +
                std::to_string(mie_phase_function_g) + ",\n" +
                density_profile(absorption_density) + ",\n" +
                to_string(
                    absorption_extinction, lambdas, length_unit_in_meters) + ",\n" +
                to_string(ground_albedo, lambdas, 1.0) + ",\n" +
                std::to_string(cos(max_sun_zenith_angle)) + ");\n" +
                functions_glsl;
    };
    transmittance_texture_ = new OpenGLTexture(TRANSMITTANCE_TEXTURE_WIDTH,
                                               TRANSMITTANCE_TEXTURE_HEIGHT);
    scattering_texture_ = new OpenGLTexture(
                SCATTERING_TEXTURE_WIDTH,
                SCATTERING_TEXTURE_HEIGHT,
                SCATTERING_TEXTURE_DEPTH,
                combine_scattering_textures||!rgb_format_supported_?GL_RGBA:GL_RGB,
                half_precision);
    if(!combine_scattering_textures){
        optional_single_mie_scattering_texture_ =
                new OpenGLTexture(
                    SCATTERING_TEXTURE_WIDTH,
                    SCATTERING_TEXTURE_HEIGHT,
                    SCATTERING_TEXTURE_DEPTH,
                    rgb_format_supported_ ? GL_RGB : GL_RGBA,
                    half_precision);
    }
    else
        optional_single_mie_scattering_texture_ = nullptr;
    irradiance_texture_ = new OpenGLTexture(
                IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
    std::string kAtmosphereShader =
            getShaderFromFile(":/AtmosphereModel/glsl/kAtmosphereShader.frag")
            .toStdString();
    std::string shader = glsl_header_factory_({kLambdaR, kLambdaG, kLambdaB})
            + "\n" + kAtmosphereShader;

    // 编译着色器
    atmosphere_shader_ = new QOpenGLShader(QOpenGLShader::Fragment);
    atmosphere_shader_->compileSourceCode(QString::fromStdString(shader));
    QString logger;
    logger = atmosphere_shader_->log();
    if(!logger.isEmpty()){
        qDebug() << "Atmosphere Shader Error::\n"
                 << logger;
    }

    Quad m_quad;
    quad = new Entity(m_quad);
}

Atmosphere::AtmosphereModel::~AtmosphereModel(){
    if(quad)delete quad;
    quad = nullptr;
    if(transmittance_texture_)
        delete transmittance_texture_;
    if(scattering_texture_)
        delete scattering_texture_;
    if(optional_single_mie_scattering_texture_)
        delete optional_single_mie_scattering_texture_;
    if(irradiance_texture_)
        delete irradiance_texture_;
    if(atmosphere_shader_)
        delete atmosphere_shader_;

}

void Atmosphere::AtmosphereModel::init(unsigned int num_scattering_orders){
    // 分配一些中间暂存纹理
    OpenGLTexture delta_irradiance_texture(
                IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
    OpenGLTexture delta_rayleigh_scattering_texture(
                SCATTERING_TEXTURE_WIDTH,
                SCATTERING_TEXTURE_HEIGHT,
                SCATTERING_TEXTURE_DEPTH,
                rgb_format_supported_ ? GL_RGB : GL_RGBA,
                half_precision_);
    OpenGLTexture delta_mie_scattering_texture(
                SCATTERING_TEXTURE_WIDTH,
                SCATTERING_TEXTURE_HEIGHT,
                SCATTERING_TEXTURE_DEPTH,
                rgb_format_supported_ ? GL_RGB : GL_RGBA,
                half_precision_);
    OpenGLTexture delta_scattering_density_texture(
                SCATTERING_TEXTURE_WIDTH,
                SCATTERING_TEXTURE_HEIGHT,
                SCATTERING_TEXTURE_DEPTH,
                rgb_format_supported_ ? GL_RGB : GL_RGBA,
                half_precision_);
    // delta_multiple_scattering_texture用于计算3次及3次以上的多重散射
    // 而delta_rayleigh_scattering_texture仅用于双重散射,为了节省内存，公用同一个纹理
    OpenGLTexture delta_multiple_scattering_texture(
                delta_rayleigh_scattering_texture.getTextureId());

    // 创建帧缓冲
    F(OGL_Function);
    GLuint fbo;
    f->glGenFramebuffers(1, &fbo);
    f->glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    // 预计算到纹理中
    vec3 lambdas{kLambdaR, kLambdaG, kLambdaB};
    precompute(delta_irradiance_texture,
               delta_rayleigh_scattering_texture,
               delta_mie_scattering_texture,
               delta_scattering_density_texture,
               delta_multiple_scattering_texture,
               lambdas,
               num_scattering_orders);

    // 释放一些资源
    f->glUseProgram(0);
    f->glBindFramebuffer(GL_FRAMEBUFFER, 0);
    f->glDeleteFramebuffers(1, &fbo);
    assert(f->glGetError() == 0);
}

void Atmosphere::AtmosphereModel::setProgramUniforms(
        GLuint program,
        GLuint transmittance_texture_unit,
        GLuint scattering_texture_unit,
        GLuint irradiance_texture_unit,
        GLuint single_mie_scattering_texture_unit) const{
    // 绑定四个纹理，用于渲染的需求
    F(OGL_Function);
    f->glActiveTexture(GL_TEXTURE0 + transmittance_texture_unit);
    f->glBindTexture(GL_TEXTURE_2D, transmittance_texture_->getTextureId());

    f->glActiveTexture(GL_TEXTURE0 + scattering_texture_unit);
    f->glBindTexture(GL_TEXTURE_3D, scattering_texture_->getTextureId());

    f->glActiveTexture(GL_TEXTURE0 + irradiance_texture_unit);
    f->glBindTexture(GL_TEXTURE_2D, irradiance_texture_->getTextureId());

    int location = f->glGetUniformLocation(program, "transmittance_texture");
    if(location == -1)qDebug() << "transmittance_texture error!";
    f->glUniform1i(location,transmittance_texture_unit);

    location = f->glGetUniformLocation(program, "scattering_texture");
    if(location == -1)qDebug() << "scattering_texture error!";
    f->glUniform1i(location,scattering_texture_unit);

    location = f->glGetUniformLocation(program, "irradiance_texture");
    if(location == -1)qDebug() << "irradiance_texture error!";
    f->glUniform1i(location,irradiance_texture_unit);

    if (optional_single_mie_scattering_texture_ != nullptr) {
        f->glActiveTexture(GL_TEXTURE0 + single_mie_scattering_texture_unit);
        f->glBindTexture(GL_TEXTURE_3D, optional_single_mie_scattering_texture_->getTextureId());
        location = f->glGetUniformLocation(program, "single_mie_scattering_texture");
        if(location == -1)qDebug() << "single_mie_scattering_texture error!";
        f->glUniform1i(location,single_mie_scattering_texture_unit);
    }
}

void Atmosphere::AtmosphereModel::precompute(
        OpenGLTexture &delta_irradiance_texture,
        OpenGLTexture &delta_rayleigh_scattering_texture,
        OpenGLTexture &delta_mie_scattering_texture,
        OpenGLTexture &delta_scattering_density_texture,
        OpenGLTexture &delta_multiple_scattering_texture,
        const Atmosphere::AtmosphereModel::vec3 &lambdas,
        unsigned int num_scattering_orders) {
    // 加载着色器代码
    std::string header = glsl_header_factory_(lambdas);
    QString directory(":/AtmosphereModel/glsl/");
    std::string kVertexShader =
            getShaderFromFile(QString(directory+"kVertexShader.vert")).toStdString();
    std::string kComputeTransmittanceShader =
            getShaderFromFile(directory+"kComputeTransmittanceShader.frag").toStdString();
    std::string kComputeDirectIrradianceShader =
            getShaderFromFile(directory+"kComputeDirectIrradianceShader.frag").toStdString();
    std::string kGeometryShader =
            getShaderFromFile(directory+"kGeometryShader.geom").toStdString();
    std::string kComputeSingleScatteringShader =
            getShaderFromFile(directory+"kComputeSingleScatteringShader.frag").toStdString();
    std::string kComputeScatteringDensityShader =
            getShaderFromFile(directory+"kComputeScatteringDensityShader.frag").toStdString();
    std::string kComputeIndirectIrradianceShader =
            getShaderFromFile(directory+"kComputeIndirectIrradianceShader.frag").toStdString();
    std::string kComputeMultipleScatteringShader =
            getShaderFromFile(directory+"kComputeMultipleScatteringShader.frag").toStdString();

    // 编译着色器代码
    ShaderProgram compute_transmittance(
                kVertexShader.c_str(),
                (header+kComputeTransmittanceShader).c_str());
    ShaderProgram compute_direct_irradiance(
                kVertexShader.c_str(),
                (header+kComputeDirectIrradianceShader).c_str());
    ShaderProgram compute_single_scattering(
                kVertexShader.c_str(),
                (header+kComputeSingleScatteringShader).c_str(),
                kGeometryShader.c_str()
                );
    ShaderProgram compute_scattering_density(
                kVertexShader.c_str(),
                (header+kComputeScatteringDensityShader).c_str(),
                kGeometryShader.c_str());
    ShaderProgram compute_indirect_irradiance(
                kVertexShader.c_str(),
                (header+kComputeIndirectIrradianceShader).c_str());
    ShaderProgram compute_multiple_scattering(
                kVertexShader.c_str(),
                (header+kComputeMultipleScatteringShader).c_str(),
                kGeometryShader.c_str());

    const GLuint kDrawBuffers[4] = {
        GL_COLOR_ATTACHMENT0,
        GL_COLOR_ATTACHMENT1,
        GL_COLOR_ATTACHMENT2,
        GL_COLOR_ATTACHMENT3
    };

    // 开始预计算
    F(OGL_Function);
    // 设置RBG和Alpha的混合模式,用于实现叠加操作
    f->glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    f->glBlendFuncSeparate(GL_ONE, GL_ONE, GL_ONE, GL_ONE);

    // 计算光学长度，并将其存储到transmittance_texture_中
    {
        f->glFramebufferTexture(
                    GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                    transmittance_texture_->getTextureId(), 0);
        f->glDrawBuffer(GL_COLOR_ATTACHMENT0);
        f->glViewport(0, 0, TRANSMITTANCE_TEXTURE_WIDTH,
                      TRANSMITTANCE_TEXTURE_HEIGHT);
        compute_transmittance.use();
        drawQuad({});
        compute_transmittance.release();
    }

    // 计算直接辐照度，存储到delta_irradiance_texture中
    // irradiance_texture_存储地面接收的天空辐照度
    {
        f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                delta_irradiance_texture.getTextureId(), 0);
        //f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1,
        //                        irradiance_texture_->getTextureId(), 0);
        f->glDrawBuffers(1, kDrawBuffers);
        f->glViewport(0, 0, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
        compute_direct_irradiance.use();
        transmittance_texture_->bind(0);
        compute_direct_irradiance.setInt("transmittance_texture",0);
        drawQuad({false});
        compute_direct_irradiance.release();
    }

    // 计算rayleigh单次内散射积分和mie单次内散射积分
    // 存储到delta_rayleigh_scattering_texture和delta_mie_scattering_texture中
    {
        f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                delta_rayleigh_scattering_texture.getTextureId(), 0);
        f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1,
                                delta_mie_scattering_texture.getTextureId(), 0);
        f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2,
                                scattering_texture_->getTextureId(), 0);
        // 是否将mie和rayleigh计算结果合并到一个纹理
        if (optional_single_mie_scattering_texture_ != nullptr) {
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3,
                                    optional_single_mie_scattering_texture_->getTextureId(),
                                    0);
            f->glDrawBuffers(4, kDrawBuffers);
        }
        else {
            f->glDrawBuffers(3, kDrawBuffers);
        }
        f->glViewport(0, 0, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
        compute_single_scattering.use();
        transmittance_texture_->bind(0);
        compute_single_scattering.setInt("transmittance_texture",0);
        for (unsigned int layer = 0; layer < SCATTERING_TEXTURE_DEPTH; ++layer) {
            compute_single_scattering.setInt("layer",layer);// 一层一层计算到3D纹理中
            drawQuad({false, false, false, false});
        }
        compute_single_scattering.release();
    }

    // 依次计算第2重、第3重...第n重散射
    {
        for (unsigned int scattering_order = 2;
             scattering_order <= num_scattering_orders; ++scattering_order) {
            // 计算特定(r,mu)接收的辐照度，存储到delta_scattering_density_texture中
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                    delta_scattering_density_texture.getTextureId(), 0);
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, 0, 0);
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, 0, 0);
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, 0, 0);
            f->glDrawBuffer(GL_COLOR_ATTACHMENT0);
            f->glViewport(0, 0, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
            compute_scattering_density.use();
            transmittance_texture_->bind(0);
            delta_rayleigh_scattering_texture.bind3D(1);
            delta_mie_scattering_texture.bind3D(2);
            delta_multiple_scattering_texture.bind3D(3);
            delta_irradiance_texture.bind(4);
            compute_scattering_density.setInt("transmittance_texture",0);
            compute_scattering_density.setInt("single_rayleigh_scattering_texture",1);
            compute_scattering_density.setInt("single_mie_scattering_texture",2);
            compute_scattering_density.setInt("multiple_scattering_texture",3);
            compute_scattering_density.setInt("irradiance_texture",4);
            compute_scattering_density.setInt("scattering_order",scattering_order);
            for (unsigned int layer = 0; layer < SCATTERING_TEXTURE_DEPTH; ++layer) {
                compute_scattering_density.setInt("layer", layer);
                drawQuad({});
            }
            compute_scattering_density.release();

            // 计算间接辐照度，存储到delta_irradiance_texture中,
            // 然后累加到irradiance_texture_
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                    delta_irradiance_texture.getTextureId(), 0);
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1,
                                    irradiance_texture_->getTextureId(), 0);
            f->glDrawBuffers(2, kDrawBuffers);
            f->glViewport(0, 0, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
            compute_indirect_irradiance.use();
            compute_indirect_irradiance.setInt("single_rayleigh_scattering_texture",0);
            compute_indirect_irradiance.setInt("single_mie_scattering_texture",1);
            compute_indirect_irradiance.setInt("multiple_scattering_texture",2);
            compute_indirect_irradiance.setInt("scattering_order",scattering_order - 1);
            delta_rayleigh_scattering_texture.bind3D(0);
            delta_mie_scattering_texture.bind3D(1);
            delta_multiple_scattering_texture.bind3D(2);
            drawQuad({false, true});
            compute_indirect_irradiance.release();

            // 计算多重散射，存储到delta_multiple_scattering_texture
            // 累加到scattering_texture_
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                    delta_multiple_scattering_texture.getTextureId(), 0);
            f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1,
                                    scattering_texture_->getTextureId(), 0);
            f->glDrawBuffers(2, kDrawBuffers);
            f->glViewport(0, 0, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
            compute_multiple_scattering.use();
            transmittance_texture_->bind(0);
            delta_scattering_density_texture.bind3D(1);
            compute_multiple_scattering.setInt("transmittance_texture",0);
            compute_multiple_scattering.setInt("scattering_density_texture",1);
            for (unsigned int layer = 0; layer < SCATTERING_TEXTURE_DEPTH; ++layer) {
                compute_multiple_scattering.setInt("layer", layer);
                drawQuad({false, true});
            }
            compute_multiple_scattering.release();
        }
    }

    // 将颜色缓冲附件都初始为原样
    f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, 0, 0);
    f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, 0, 0);
    f->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, 0, 0);
}

void Atmosphere::AtmosphereModel::drawQuad(const std::vector<bool> &enable_blend){
    // 绘制，仅仅是为了precomputation
    F(OGL_Function);
    for (unsigned int i = 0; i < enable_blend.size(); ++i) {
        if (enable_blend[i]) {
            f->glEnablei(GL_BLEND, i);
        }
    }
    quad->draw(nullptr);
    for (unsigned int i = 0; i < enable_blend.size(); ++i) {
        f->glDisablei(GL_BLEND, i);
    }
}

double Atmosphere::AtmosphereModel::Interpolate(
        const std::vector<double> &wavelengths,
        const std::vector<double> &wavelength_function, double wavelength){
    if (wavelength < wavelengths[0]) {
        return wavelength_function[0];
    }
    for (unsigned int i = 0; i < wavelengths.size() - 1; ++i) {
        if (wavelength < wavelengths[i + 1]) {// 线性插值
            double u =
                    (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
            return
                    wavelength_function[i] * (1.0 - u) + wavelength_function[i + 1] * u;
        }
    }
    return wavelength_function[wavelength_function.size() - 1];
}

bool Atmosphere::AtmosphereModel::IsFramebufferRgbFormatSupported(
        bool half_precision){
    // 检测是否支持RGB格式
    F(OGL_Function);
    GLuint test_fbo = 0;
    f->glGenFramebuffers(1, &test_fbo);
    f->glBindFramebuffer(GL_FRAMEBUFFER, test_fbo);
    GLuint test_texture = 0;
    f->glGenTextures(1, &test_texture);
    f->glBindTexture(GL_TEXTURE_2D, test_texture);
    f->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    f->glTexImage2D(GL_TEXTURE_2D, 0, half_precision ? GL_RGB16F : GL_RGB32F,
                    1, 1, 0, GL_RGB, GL_FLOAT, NULL);
    f->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                              GL_TEXTURE_2D, test_texture, 0);
    bool rgb_format_supported =
            f->glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    f->glDeleteTextures(1, &test_texture);
    f->glDeleteFramebuffers(1, &test_fbo);
    return rgb_format_supported;
}

QString Atmosphere::AtmosphereModel::getShaderFromFile(QString path){
    // 从文件读取shader代码
    QFile file(path);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QByteArray con = file.readAll();
    QString context(con);
    file.close();
    return context;
}


