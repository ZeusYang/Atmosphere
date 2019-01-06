#ifndef THIRDPERSONCAMERA_H
#define THIRDPERSONCAMERA_H

#include "Utilities/Camera3D.h"
#include "Utilities/Player.h"

class ThirdPersonCamera : public Camera3D{
public:
    ThirdPersonCamera(Player* p);
    virtual ~ThirdPersonCamera();

    void move();
    virtual QVector3D getPosition();

private:
    float yaw,pitch;
    float angleAroundPlayer,distanceFromPlayer;
    QVector3D position;
    Player* player;

    virtual const QMatrix4x4& toMatrix();

    void calculateZoom();
    void calculatePitch();
    void calculateAngleAroundPlayer();
    float calculateHorizontalDistance();
    float calculateVerticalDistance();
    void calculateCameraPosition(float horizDistance, float verticDistance);
};

inline QVector3D ThirdPersonCamera::getPosition(){
    return this->position;
}

#endif // THIRDPERSONCAMERA_H
