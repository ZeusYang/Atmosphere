#include "GodCamera.h"
#include "Utilities/InputManager.h"
#include <cmath>

GodCamera::GodCamera(double z, double a, double d)
    :view_zenith_angle_radians_(z),
      view_azimuth_angle_radians_(a),
      view_distance_meters_(d)
{
    update();
}

GodCamera::~GodCamera() {}

void GodCamera::setProjectMatrix(float fov, float asp, float zNea, float zFar){
    m_invProjectMatrix.setToIdentity();
    m_invProjectMatrix.perspective(fov, asp, zNea, zFar);
    m_invProjectMatrix = m_invProjectMatrix.inverted();
}

void GodCamera::updateViewMatrix(double dz, double da, double dd){
    view_zenith_angle_radians_  = dz;
    view_azimuth_angle_radians_ = da;
    view_distance_meters_       = dd;
    update();
}

void GodCamera::update(){
    float cos_z = std::cos(view_zenith_angle_radians_);
    float sin_z = std::sin(view_zenith_angle_radians_);
    float cos_a = std::cos(view_azimuth_angle_radians_);
    float sin_a = std::sin(view_azimuth_angle_radians_);
    float ux[3] = { -sin_a, cos_a, 0.0 };
    float uy[3] = { -cos_z * cos_a, -cos_z * sin_a, sin_z };
    float uz[3] = { sin_z * cos_a, sin_z * sin_a, cos_z };
    float l     = view_distance_meters_;

    float model_from_view[16] = {
        ux[0], uy[0], uz[0], uz[0] * l,
        ux[1], uy[1], uz[1], uz[1] * l,
        ux[2], uy[2], uz[2], uz[2] * l,
        0.0, 0.0, 0.0, 1.0
    };

    position = QVector3D(model_from_view[3],model_from_view[7],model_from_view[11]);

    m_invViewMatrix = QMatrix4x4(model_from_view);
    m_invViewMatrix = m_invViewMatrix.transposed();
}
