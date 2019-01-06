#include "Player.h"
#include <QtMath>

Player::Player(){
    rotY = rotX = rotZ = 0;
    position = QVector3D(0,0,0);
    matrix.setTranslation(position);

    currentSpeed = 0;
    currentTurnSpeed = 0;
    upwardsSpeed = 0;
    isInAir = false;
    m_frameTimer.start();
}

Player::Player(Entity *entity, QVector3D pos) : target(entity){
    rotY = rotX = rotZ = 0;
    position = pos;
    matrix.setTranslation(position);

    currentSpeed = 0;
    currentTurnSpeed = 0;
    upwardsSpeed = 0;
    isInAir = false;
    m_frameTimer.start();
}

Player::Player(Player &target){
    rotX = target.rotX;
    rotY = target.rotY;
    rotZ = target.rotZ;
    position = target.getPosition();
    matrix.setTranslation(position);

    currentSpeed = 0;
    currentTurnSpeed = 0;
    upwardsSpeed = 0;
    isInAir = false;
    m_frameTimer.start();
}

Player &Player::operator=(Player &target){
    rotX = target.rotX;
    rotY = target.rotY;
    rotZ = target.rotZ;
    position = target.getPosition();
    matrix.setTranslation(position);

    currentSpeed = 0;
    currentTurnSpeed = 0;
    upwardsSpeed = 0;
    isInAir = false;
    m_frameTimer.start();
}

void Player::move(){
    // calculate delta time
    qint64 elapsedTime = m_frameTimer.elapsed();
    float deltaTime = (elapsedTime)/1000.0f;

    // check key board operation
    checkInputs();

    // rotate the angle with y-axis.
    matrix.rotate(currentTurnSpeed*deltaTime,0.0,1.0,0.0);
    rotY += currentTurnSpeed*deltaTime;
    rotY = fmod(rotY,360.0f);

    // calc the delta distance
    float distance = currentSpeed*deltaTime;
    float dx = (float) (distance * qSin(qDegreesToRadians(rotY)));
    float dz = (float) (distance * qCos(qDegreesToRadians(rotY)));

    // jump it
    upwardsSpeed += GRAVITY*deltaTime;
    position += QVector3D(dx,upwardsSpeed,dz);
    if(position.y() < 0){
        upwardsSpeed = 0;
        isInAir = false;
        position.setY(0);
    }
    matrix.setTranslation(position);

    // Set entity model matrix
    if(target){
        QMatrix4x4 record = matrix.toMatrix();
        target->setModelMatrix(record);
    }
    m_frameTimer.start();
}

QMatrix4x4 Player::toMatrix(){
    return matrix.toMatrix();
}

QVector3D Player::getPosition() const{
    // Player's position
    return this->position;
}

float Player::getRotY() const{
    return this->rotY;
}

void Player::checkInputs(){
    if(InputManager::keyPressed(Qt::Key_W)){
        currentSpeed = RUN_SPEED;
    }
    else if(InputManager::keyPressed(Qt::Key_S)){
        currentSpeed = -RUN_SPEED;
    }
    else{
        currentSpeed = 0;
    }

    if(InputManager::keyPressed(Qt::Key_D)){
        currentTurnSpeed = -TURN_SPEED;
    }
    else if(InputManager::keyPressed(Qt::Key_A)){
        currentTurnSpeed = TURN_SPEED;
    }
    else{
        currentTurnSpeed = 0;
    }

    if(InputManager::keyPressed(Qt::Key_Space)){
        jump();
    }
}

void Player::jump() {
    if (!isInAir) {
        upwardsSpeed = JUMP_POWER;
        isInAir = true;
    }
}
