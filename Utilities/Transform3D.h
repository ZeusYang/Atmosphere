#ifndef TRANSFORM3D_H
#define TRANSFORM3D_H

#include <QVector3D>
#include <QMatrix4x4>
#include <QQuaternion>

class Transform3D {
public:
    static const QVector3D LocalForward;
    static const QVector3D LocalUp;
    static const QVector3D LocalRight;

    // Constructors
    Transform3D();

    // Transform By (Add/Scale)
    void translate(const QVector3D &dt);
    void translate(float dx, float dy, float dz);
    void scale(const QVector3D &ds);
    void scale(float dx, float dy, float dz);
    void scale(float factor);
    void rotate(const QQuaternion &dr);
    void rotate(float angle, const QVector3D &axis);
    void rotate(float angle, float ax, float ay, float az);

    // Transform To (Setters)
    void setTranslation(const QVector3D &t);
    void setTranslation(float x, float y, float z);
    void setScale(const QVector3D &s);
    void setScale(float x, float y, float z);
    void setScale(float k);
    void setRotation(const QQuaternion &r);
    void setRotation(float angle, const QVector3D &axis);
    void setRotation(float angle, float ax, float ay, float az);

    // Accessors
    const QVector3D& scale() const;
    const QVector3D& translation() const;
    const QQuaternion& rotation() const;
    QMatrix4x4& toMatrix();

    // Queries
    QVector3D up() const;
    QVector3D forward() const;
    QVector3D right() const;

private:
    bool m_dirty;
    QVector3D m_scale;
    QVector3D m_translation;
    QQuaternion m_rotation;
    QMatrix4x4 m_world;
};

inline Transform3D::Transform3D() : m_dirty(true), m_scale(1.0f, 1.0f, 1.0f) {}

// Transform By (Add/Scale)
inline void Transform3D::translate(float dx, float dy,float dz) { translate(QVector3D(dx, dy, dz)); }
inline void Transform3D::scale(float dx, float dy,float dz) { scale(QVector3D(dx, dy, dz)); }
inline void Transform3D::scale(float factor) { scale(QVector3D(factor, factor, factor)); }
inline void Transform3D::rotate(float angle, const QVector3D &axis) { rotate(QQuaternion::fromAxisAndAngle(axis, angle)); }
inline void Transform3D::rotate(float angle, float ax, float ay,float az) { rotate(QQuaternion::fromAxisAndAngle(ax, ay, az, angle)); }

// Transform To (Setters)
inline void Transform3D::setTranslation(float x, float y, float z) { setTranslation(QVector3D(x, y, z)); }
inline void Transform3D::setScale(float x, float y, float z) { setScale(QVector3D(x, y, z)); }
inline void Transform3D::setScale(float k) { setScale(QVector3D(k, k, k)); }
inline void Transform3D::setRotation(float angle, const QVector3D &axis) { setRotation(QQuaternion::fromAxisAndAngle(axis, angle)); }
inline void Transform3D::setRotation(float angle, float ax, float ay, float az) { setRotation(QQuaternion::fromAxisAndAngle(ax, ay, az, angle)); }

// Accessors
inline const QVector3D& Transform3D::translation() const { return m_translation; }
inline const QVector3D& Transform3D::scale() const { return m_scale; }
inline const QQuaternion& Transform3D::rotation() const { return m_rotation; }

#endif // TRANSFORM3D_H
