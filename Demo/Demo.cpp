#include "Demo.h"
#include "Utilities/Geometry/Quad.h"
#include <QOpenGLShader>
#include <QString>
#include <assert.h>
#include <cmath>
#include "Utilities/OpenGLVersion.h"
#include "Utilities/InputManager.h"

using namespace Atmosphere;


#define F(c) c *f = (GlobalContext::contextFunc)

constexpr double kPi = 3.1415926;
// 太阳角半径
constexpr double kSunAngularRadius = 0.00935 / 2.0;
// 太阳立体角
constexpr double kSunSolidAngle = kPi * kSunAngularRadius * kSunAngularRadius;
// 长度单位
constexpr double kLengthUnitInMeters = 1000.0;

Demo::Demo():
    use_constant_solar_spectrum_(false),
    use_ozone_(true),
    use_combined_textures_(true),
    use_half_precision_(true),
    view_distance_meters_(9000.0),
    view_zenith_angle_radians_(1.47),
    view_azimuth_angle_radians_(-0.1),
    sun_zenith_angle_radians_(1.3),
    sun_azimuth_angle_radians_(2.9),
    exposure_(10.0),
    _program(nullptr),
    model(nullptr),
    light_shaft(1),
    num_scattering_orders(5),
    use_rayleigh_scattering(true),
    use_mie_scattering(true),
    mie_g_coefficient(0.8)
{
    // 创建一个屏幕大小的四边形
    Quad m_quad;
    screen = new Entity(m_quad);
    // 创建大气层模型
    initModel();
    // 摄像机
    camera = new GodCamera(view_zenith_angle_radians_,
                           view_azimuth_angle_radians_,
                           view_distance_meters_/kLengthUnitInMeters);
}

Demo::~Demo(){
    if(_program)delete _program;
    if(model)delete model;
    if(screen)delete screen;
    if(camera)delete camera;
    _program = nullptr;
    model  = nullptr;
    screen = nullptr;
    camera = nullptr;
}

void Demo::draw(){
    if(!_program)return;
    // 渲染场景
    _program->use();
    _program->setVector3("camera",camera->getPosition());
    _program->setFloat("exposure",exposure_);
    _program->setInt("light_shaft",light_shaft);
    _program->setMatrix4("model_from_view",camera->invViewMatrix());
    _program->setMatrix4("view_from_clip",camera->invProjectMatrix());
    _program->setVector3("sun_direction",QVector3D(
                             cos(sun_azimuth_angle_radians_) * sin(sun_zenith_angle_radians_),
                             sin(sun_azimuth_angle_radians_) * sin(sun_zenith_angle_radians_),
                             cos(sun_zenith_angle_radians_)));
    screen->draw(nullptr);
    _program->release();
}

void Demo::resize(float width, float height){
    camera->setProjectMatrix(50.0f,(float)width/(float)height,0.5f,1000.0f);
}

void Demo::update(){
    constexpr double kScale = 500.0;
    // 移动太阳
    if(InputManager::keyPressed(Qt::Key_Control) &&
            InputManager::buttonPressed(Qt::LeftButton)){
        sun_zenith_angle_radians_ +=
                (double)InputManager::mouseDelta().y()*0.5f/kScale;
        sun_zenith_angle_radians_ =
                std::max(0.0, std::min(kPi, sun_zenith_angle_radians_));
        sun_azimuth_angle_radians_ -=
                (double)InputManager::mouseDelta().x()*0.5f/kScale;
    }
    // 变换视角
    else if(InputManager::buttonPressed(Qt::LeftButton)){
        view_zenith_angle_radians_ -=
                (double)InputManager::mouseDelta().y()*0.5f/kScale;
        view_zenith_angle_radians_ =
                std::max(0.0, std::min(kPi / 2.0, view_zenith_angle_radians_));
        view_azimuth_angle_radians_ -=
                (double)InputManager::mouseDelta().x()*0.5f/kScale;
    }
    // 移动距离
    float wheelDelta = InputManager::wheelDelta().y();
    if(wheelDelta < 0){
        view_distance_meters_ *= 1.05;
    }
    else if(wheelDelta > 0){
        view_distance_meters_ /= 1.05;
    }
    // 更新视角矩阵
    camera->updateViewMatrix(view_zenith_angle_radians_,
                             view_azimuth_angle_radians_,
                             view_distance_meters_/kLengthUnitInMeters);
    InputManager::setWheelDelta(QPoint(0,0));
}

void Demo::setExposure(double expos){
    exposure_ = expos;
}

bool Demo::setSpectrum(bool use){
    if(use_constant_solar_spectrum_ == use)return false;
    use_constant_solar_spectrum_ = use;
    return true;
}

bool Demo::setOZone(bool use){
    if(use_ozone_ == use)return false;
    use_ozone_ = use;
    return true;
}

bool Demo::setCombined(bool use){
    if(use_combined_textures_ == use)return false;
    use_combined_textures_ = use;
    return true;
}

bool Demo::setHalfpre(bool use){
    if(use_half_precision_ == use)return false;
    use_half_precision_ = use;
    return true;
}

void Demo::setLightshaft(bool use){
    if(use)light_shaft = 1;
    else light_shaft = 0;
}

bool Demo::setNumScattering(int nums){
    if(num_scattering_orders == nums)return false;
    num_scattering_orders = nums;
    return true;
}

bool Demo::setRayleigh(bool use){
    if(use == use_rayleigh_scattering)return false;
    use_rayleigh_scattering = use;
    return true;
}

bool Demo::setMie(bool use){
    if(use == use_mie_scattering)return false;
    use_mie_scattering = use;
    return true;
}

void Demo::initModel(){
    if(_program)delete _program;
    if(model)delete model;
    _program = nullptr;
    model = nullptr;
    // 太阳光谱
    constexpr int kLambdaMin = 360;
    constexpr int kLambdaMax = 830;
    constexpr double kSolarIrradiance[48] = {
        1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.68870, 1.61253,
        1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.82980,
        1.86850, 1.89310, 1.85149, 1.85040, 1.83410, 1.83450, 1.81470, 1.78158, 1.7533,
        1.69650, 1.68194, 1.64654, 1.60480, 1.52143, 1.55622, 1.51130, 1.47400, 1.4482,
        1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.23670, 1.20820,
        1.18737, 1.14683, 1.12362, 1.10580, 1.07124, 1.04992
    };
    // http://www.iup.uni-bremen.de/gruppen/molspec/databases
    // /referencespectra/o3spectra2011/index.html
    constexpr double kOzoneCrossSection[48] = {
        1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
        8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
        1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
        4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
        2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
        6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
        2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
    };
    // https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
    constexpr double kDobsonUnit = 2.687e20;
    // 臭氧层最大密度
    constexpr double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;
    // 与波长无关的太阳辐照度光谱
    constexpr double kConstantSolarIrradiance = 1.5;
    constexpr double kBottomRadius = 6360000.0;
    constexpr double kTopRadius = 6420000.0;
    constexpr double kRayleigh = 1.24062e-6;
    // rayleigh散射缩放高度
    constexpr double kRayleighScaleHeight = 8000.0;
    // mie散射缩放高度
    constexpr double kMieScaleHeight = 1200.0;
    constexpr double kMieAngstromAlpha = 0.0;
    constexpr double kMieAngstromBeta = 5.328e-3;
    constexpr double kMieSingleScatteringAlbedo = 0.9;
    double kMiePhaseFunctionG = mie_g_coefficient;
    constexpr double kGroundAlbedo = 0.1;
    const double max_sun_zenith_angle =
            (use_half_precision_ ? 102.0 : 120.0) / 180.0 * kPi;

    // rayleigh层，即空气分子层     //width,exp_term,exp_scale,linear_term,constant_term
    DensityProfileLayer rayleigh_layer(0.0, 1.0*use_rayleigh_scattering,
                                       -1.0/kRayleighScaleHeight*use_rayleigh_scattering,
                                       0.0, 0.0);
    // mie层，气溶胶层            //width,exp_term,exp_scale,linear_term,constant_term
    DensityProfileLayer mie_layer(0.0, 1.0*use_mie_scattering,
                                  -1.0/kMieScaleHeight*use_mie_scattering, 0.0, 0.0);
    std::vector<DensityProfileLayer> ozone_density;
    ozone_density.push_back(
                DensityProfileLayer(25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0));
    ozone_density.push_back(
                DensityProfileLayer(0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0));

    std::vector<double> wavelengths;            // 波长
    std::vector<double> solar_irradiance;       // 太阳辐照度
    std::vector<double> rayleigh_scattering;    // rayleigh散射
    std::vector<double> mie_scattering;         // mie散射
    std::vector<double> mie_extinction;         // mie消光
    std::vector<double> absorption_extinction;  // 吸收光线的空气分子消光
    std::vector<double> ground_albedo;          // 地面反照率
    for (int l = kLambdaMin; l <= kLambdaMax; l += 10) {
        double lambda = static_cast<double>(l) * 1e-3;
        double mie =
                kMieAngstromBeta / kMieScaleHeight * pow(lambda, -kMieAngstromAlpha);
        // 太阳光波波长
        wavelengths.push_back(l);
        // 太阳辐照度
        if (use_constant_solar_spectrum_) {//常量
            solar_irradiance.push_back(kConstantSolarIrradiance);
        } else {
            solar_irradiance.push_back(kSolarIrradiance[(l - kLambdaMin) / 10]);
        }
        rayleigh_scattering.push_back(kRayleigh * pow(lambda, -4));
        mie_scattering.push_back(mie * kMieSingleScatteringAlbedo);
        mie_extinction.push_back(mie);
        absorption_extinction.push_back(
                    use_ozone_*kMaxOzoneNumberDensity * kOzoneCrossSection[(l - kLambdaMin) / 10]);
        ground_albedo.push_back(kGroundAlbedo);
    }
    // 创建新模型
    model = new AtmosphereModel(
                wavelengths,                            // 太阳波长，单位nm
                solar_irradiance,                       // 太阳辐照度
                kSunAngularRadius,                      // 太阳角半径
                kBottomRadius,                          // 大气层底层到星球中心的距离(内半径)
                kTopRadius,                             // 大气层外层到星球中心的距离(外半径)
                {rayleigh_layer},                       // 大气空气分子密度分布
                rayleigh_scattering,                    // rayleigh散射系数=海拔高度h处的rayleigh_scattering*rayleigh_density
                {mie_layer},                            // 大气气溶胶密度分布
                mie_scattering,                         // (气溶胶)mie散射系数=海拔高度h处的mie_scattering*mie_density
                mie_extinction,                         // (气溶胶)mie消光系数=海拔高度h处mie_extinction*mie_density
                kMiePhaseFunctionG,                     // 气溶胶Cornette-Shanks的相位函数参数值g
                ozone_density,                          // 大气中吸收光线的空气分子密度
                absorption_extinction,                  // 大气中吸收光线的消光系数=海拔高度h处absorption_extinction*absorption_density
                ground_albedo,                          // 地面的平均反照率
                max_sun_zenith_angle,                   // 太阳最大的天顶角,弧度制
                kLengthUnitInMeters,                    // 长度单位
                use_combined_textures_,                 // 是否把单次mie散射和rayleigh散射以及多次散射合并到一个纹理
                use_half_precision_);                   // 采用半精度浮点数还是单精度浮点数

    model->init(num_scattering_orders);

    // 编译demo的shaders
    std::string demo_glsl =
            ShaderProgram::getShaderFromFile(
                QString::fromStdString(":/Shaders/Demo.frag")).toStdString();
    std::string fragment_shader_str =
            "#version 330\n" +
            std::string("") +
            "const float kLengthUnitInMeters = " +
            std::to_string(kLengthUnitInMeters) + ";\n" +
            demo_glsl;
    QString logger1,logger2;
    QOpenGLShader vertex(QOpenGLShader::Vertex);
    QOpenGLShader fragment(QOpenGLShader::Fragment);
    vertex.compileSourceFile(":/Shaders/Demo.vert");
    fragment.compileSourceCode(QString::fromStdString(fragment_shader_str));
    logger1 = vertex.log();
    logger2 = fragment.log();
    if(!logger1.isEmpty())qDebug() << "Demo Vertex Shader::\n" << logger1;
    if(!logger2.isEmpty())qDebug() << "Demo Fragment Shader::\n" << logger2;

    _program = new RawShaderProgram({vertex.shaderId(),
                                     fragment.shaderId(),
                                     model->getShader()->shaderId()});

    // 设置一些uniform值
    _program->use();
    model->setProgramUniforms(_program->programId(),0,1,2,3);
    _program->setVector3("earth_center",QVector3D(0.0,
                                                  0.0,
                                                  -kBottomRadius / kLengthUnitInMeters));
    _program->setVector2("sun_size",QVector2D(std::tan(kSunAngularRadius),
                                              std::cos(kSunAngularRadius)));
    _program->release();
    F(OGL_Function);
    assert(f->glGetError() == 0);
}

bool Demo::setMieGCof(double target){
    if(mie_g_coefficient == target)return false;
    mie_g_coefficient = target;
    return true;
}

void Demo::SetView(double view_distance_meters,
                   double view_zenith_angle_radians,
                   double view_azimuth_angle_radians,
                   double sun_zenith_angle_radians,
                   double sun_azimuth_angle_radians,
                   double exposure) {
    view_distance_meters_       = view_distance_meters;
    view_zenith_angle_radians_  = view_zenith_angle_radians;
    view_azimuth_angle_radians_ = view_azimuth_angle_radians;
    sun_zenith_angle_radians_   = sun_zenith_angle_radians;
    sun_azimuth_angle_radians_  = sun_azimuth_angle_radians;
    exposure_                   = exposure;
}

