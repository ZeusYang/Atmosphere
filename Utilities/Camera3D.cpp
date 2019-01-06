#include "Camera3D.h"

const QVector3D Camera3D::LocalForward(0.0f, 0.0f, -1.0f);
const QVector3D Camera3D::LocalUp(0.0f, 1.0f, 0.0f);
const QVector3D Camera3D::LocalRight(1.0f, 0.0f, 0.0f);

Camera3D::~Camera3D() {}

QVector3D Camera3D::getPosition(){
    // Get camera position
    QVector4D tmp(0,0,0,1);
    tmp = this->m_translation * tmp;
    return QVector3D(tmp.x(),tmp.y(),tmp.z());
}

QMatrix4x4 Camera3D::getInvViewMatrix(){
    if(m_dirty){
        viewMatrix = this->toMatrix();
        viewProjectMatrix = projectMatrix * viewMatrix;
        invViewMatrix = viewMatrix.inverted();
    }
    return invViewMatrix;
}

QMatrix4x4 Camera3D::getInvProjectMatrix(){
    return invProjectMatrix;
}

void Camera3D::move(){
    //nothing
}

void Camera3D::moveTo(QVector3D position){
    // Move camera to position
    this->translate(position);
}

void Camera3D::initPerspectProject(float fov, float asp, float zNea, float zFa){
    this->fovy = fov;
    this->aspect = asp;
    this->zNear = zNea;
    this->zFar = zFa;
    projectMatrix.setToIdentity();
    projectMatrix.perspective(fovy, aspect, zNear, zFar);
    invProjectMatrix = projectMatrix.inverted();
}

void Camera3D::initOrthoProject(float left, float right, float bottom, float top, float near, float far){
    projectMatrix.setToIdentity();
    projectMatrix.ortho(left,right,bottom,top,near,far);
    invProjectMatrix = projectMatrix.inverted();
}

void Camera3D::translate(const QVector3D &dt){
    m_dirty = true;
    m_translation += dt;
}

void Camera3D::rotate(const QQuaternion &dr){
    m_dirty = true;
    m_rotation = dr * m_rotation;
}

// Transform To (Setters)
void Camera3D::setTranslation(const QVector3D &t){
    m_dirty = true;
    m_translation = t;
}

void Camera3D::setRotation(const QQuaternion &r){
    m_dirty = true;
    m_rotation = r;
}

// Accessors
const QMatrix4x4 &Camera3D::toMatrix(){
    if (m_dirty){
        m_dirty = false;
        m_world.setToIdentity();
        m_world.rotate(m_rotation.conjugate());
        m_world.translate(-m_translation);
    }
    return m_world;
}

// Queries
QVector3D Camera3D::forward() const{
    return m_rotation.rotatedVector(LocalForward);
}

QVector3D Camera3D::up() const{
    return m_rotation.rotatedVector(LocalUp);
}

QVector3D Camera3D::right() const{
    return m_rotation.rotatedVector(LocalRight);
}
