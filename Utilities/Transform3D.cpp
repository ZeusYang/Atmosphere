#include "Transform3D.h"

const QVector3D Transform3D::LocalForward(0.0f, 0.0f, 1.0f);
const QVector3D Transform3D::LocalUp(0.0f, 1.0f, 0.0f);
const QVector3D Transform3D::LocalRight(1.0f, 0.0f, 0.0f);

// Transform By (Add/Scale)
void Transform3D::translate(const QVector3D &dt){
    m_dirty = true;
    m_translation += dt;
}

void Transform3D::scale(const QVector3D &ds){
    m_dirty = true;
    m_scale *= ds;
}

void Transform3D::rotate(const QQuaternion &dr){
    m_dirty = true;
    m_rotation = dr * m_rotation;
}

// Transform To (Setters)
void Transform3D::setTranslation(const QVector3D &t){
    m_dirty = true;
    m_translation = t;
}

void Transform3D::setScale(const QVector3D &s){
    m_dirty = true;
    m_scale = s;
}

void Transform3D::setRotation(const QQuaternion &r){
    m_dirty = true;
    m_rotation = r;
}

// Accessors
QMatrix4x4 &Transform3D::toMatrix(){
    // 仅当修改了才重新计算变换矩阵
    if (m_dirty){
        m_dirty = false;
        m_world.setToIdentity();
        m_world.translate(m_translation);
        m_world.rotate(m_rotation);
        m_world.scale(m_scale);
    }
    return m_world;
}

// Queries
QVector3D Transform3D::forward() const{
    return m_rotation.rotatedVector(LocalForward);
}

QVector3D Transform3D::up() const{
    return m_rotation.rotatedVector(LocalUp);
}

QVector3D Transform3D::right() const{
    return m_rotation.rotatedVector(LocalRight);
}
