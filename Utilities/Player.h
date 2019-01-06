#ifndef PLAYER_H
#define PLAYER_H

#include "Utilities/Transform3D.h"
#include "Utilities/InputManager.h"
#include "Utilities/Entity.h"
#include <QElapsedTimer>

class Player
{
public:
    Player();
    Player(Entity *target, QVector3D pos);
    Player(Player &target);
    Player& operator=(Player &target);

    void move();
    QMatrix4x4 toMatrix();
    QVector3D getPosition()const;
    float getRotY()const;

private:
    const float RUN_SPEED = 20;
    const float TURN_SPEED = 160;
    const float GRAVITY = -20;
    const float JUMP_POWER = 5;

    float currentSpeed = 0;
    float currentTurnSpeed = 0;
    float upwardsSpeed = 0;
    bool isInAir = false;
    // 饶X、Y、Z轴旋转的角度
    float rotX,rotY,rotZ;
    Entity* target;
    QVector3D position;

    Transform3D matrix;
    QElapsedTimer m_frameTimer;

    void checkInputs();
    void jump();
};

#endif // PLAYER_H
