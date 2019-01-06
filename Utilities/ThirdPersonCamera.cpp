#include "ThirdPersonCamera.h"
#include <QtMath>

ThirdPersonCamera::ThirdPersonCamera(Player* p)
    :distanceFromPlayer(80){
    yaw = 0;
    pitch = 20;
    angleAroundPlayer = 0;
    player = p;
}

ThirdPersonCamera::~ThirdPersonCamera() {}

void ThirdPersonCamera::move(){
    calculateZoom();
    calculatePitch();
    calculateAngleAroundPlayer();
    float horizontalDistance = calculateHorizontalDistance();
    float verticalDistance = calculateVerticalDistance();
    calculateCameraPosition(horizontalDistance, verticalDistance);
    yaw = 180 - (player->getRotY() + angleAroundPlayer);
}

const QMatrix4x4& ThirdPersonCamera::toMatrix(){
    // Camera's world matrix
    m_world.setToIdentity();
    m_world.rotate(pitch,1,0,0);
    m_world.rotate(yaw,0,1,0);
    m_world.translate(-position);
    return m_world;
}

void ThirdPersonCamera::calculateCameraPosition(float horizDistance, float verticDistance){
    float theta = player->getRotY() + angleAroundPlayer;
    float offsetX = (float)(horizDistance*qSin(qDegreesToRadians(theta)));
    float offsetZ = (float)(horizDistance*qCos(qDegreesToRadians(theta)));
    position.setX(player->getPosition().x() - offsetX);
    position.setY(player->getPosition().y() + verticDistance);
    position.setZ(player->getPosition().z() - offsetZ);
}

float ThirdPersonCamera::calculateHorizontalDistance(){
    return (float)(distanceFromPlayer*qCos(qDegreesToRadians(pitch)));
}

float ThirdPersonCamera::calculateVerticalDistance(){
    return (float)(distanceFromPlayer*qSin(qDegreesToRadians(pitch)));
}

void ThirdPersonCamera::calculateZoom(){
    float zoomLevel = InputManager::wheelDelta().y() * 0.3f;
    distanceFromPlayer -= zoomLevel;

    if (distanceFromPlayer > 400)
        distanceFromPlayer = 400;
    else if (distanceFromPlayer < 1)
        distanceFromPlayer = 1;
    // set to default, otherwise it'll not stop untill to the border.
    InputManager::setWheelDelta(QPoint(0,0));
}

void ThirdPersonCamera::calculatePitch(){
    if(InputManager::buttonPressed(Qt::RightButton)){
        float pitchChange = InputManager::mouseDelta().y()*0.1f;
        pitch += pitchChange;
    }

    // limitation
    if (pitch < 1)
        pitch = 1;
    else if (pitch > 90)
        pitch = 90;
}

void ThirdPersonCamera::calculateAngleAroundPlayer(){
    if(InputManager::buttonPressed(Qt::LeftButton)){
        float angleChange = InputManager::mouseDelta().x()*0.2f;
        angleAroundPlayer -= angleChange;
    }
}
