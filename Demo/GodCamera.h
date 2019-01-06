#ifndef GODCAMERA_H
#define GODCAMERA_H

#include <QMatrix4x4>

class GodCamera {
public:
    GodCamera(double z,double a,double d);
    ~GodCamera();

    QMatrix4x4 invViewMatrix()const{return m_invViewMatrix;}
    QMatrix4x4 invProjectMatrix()const{return m_invProjectMatrix;}

    void setProjectMatrix(float fov,float asp,float zNea,float zFar);
    void updateViewMatrix(double dz,double da,double dd);

    QVector3D getPosition()const{return position;}

private:
    double view_zenith_angle_radians_;
    double view_azimuth_angle_radians_;
    double view_distance_meters_;
    QMatrix4x4 m_invViewMatrix;
    QMatrix4x4 m_invProjectMatrix;
    QVector3D position;

    void update();
};

#endif // GODCAMERA_H
