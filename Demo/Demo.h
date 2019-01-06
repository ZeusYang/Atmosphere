#ifndef DEMO_H
#define DEMO_H

#include "Utilities/Entity.h"
#include "GodCamera.h"
#include "RawShaderProgram.h"
#include "AtmosphereModel/AtmosphereModel.h"

class Demo{
public:
    Demo();
    ~Demo();

    void draw();
    void resize(float width,float height);
    void update();

    // Parameters setting
    void setExposure(double expos);
    bool setSpectrum(bool use);
    bool setOZone(bool use);
    bool setCombined(bool use);
    bool setHalfpre(bool use);
    void setLightshaft(bool use);
    bool setNumScattering(int nums);
    bool setRayleigh(bool use);
    bool setMie(bool use);
    void initModel();
    bool setMieGCof(double target);

private:
    void SetView(double view_distance_meters, double view_zenith_angle_radians,
                 double view_azimuth_angle_radians, double sun_zenith_angle_radians,
                 double sun_azimuth_angle_radians, double exposure);

    bool use_constant_solar_spectrum_;      // 设太阳光谱为常量
    bool use_ozone_;                        // 开启臭氧层
    bool use_combined_textures_;            // 开启合并纹理
    bool use_half_precision_;               // 半精度浮点数
    double exposure_;                       // 曝光度
    int  light_shaft;                       // 体积光
    unsigned int num_scattering_orders;     // 散射重数
    bool use_rayleigh_scattering;           // 开启rayleigh
    bool use_mie_scattering;                // 开启mie散射
    double mie_g_coefficient;

    double view_distance_meters_;
    double view_zenith_angle_radians_;      // 视线天顶角（弧度制)
    double view_azimuth_angle_radians_;     // 视线方位角（弧度制）
    double sun_zenith_angle_radians_;       // 太阳天顶角（弧度制）
    double sun_azimuth_angle_radians_;      // 太阳方位角（弧度制）

    Entity *screen;
    GodCamera *camera;
    RawShaderProgram *_program;
    Atmosphere::AtmosphereModel *model;
};

#endif // DEMO_H
