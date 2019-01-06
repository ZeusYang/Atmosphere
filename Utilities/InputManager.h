#ifndef INPUTMANAGER_H
#define INPUTMANAGER_H

#include <QPoint>
#include <Qt>

/**********************************
 * Input management
 **********************************/

class InputManager {
public:
    // Possible key states
    enum InputState {
        InputInvalid,
        InputRegistered,
        InputUnregistered,
        InputTriggered,
        InputPressed,
        InputReleased
    };

    // State checking
    // Keyboard input
    static InputState keyState(Qt::Key key);
    static bool keyTriggered(Qt::Key key);
    static bool keyPressed(Qt::Key key);
    static bool keyReleased(Qt::Key key);
    // Button input
    static InputState buttonState(Qt::MouseButton button);
    static bool buttonTriggered(Qt::MouseButton button);
    static bool buttonPressed(Qt::MouseButton button);
    static bool buttonReleased(Qt::MouseButton button);
    // Mouse input
    static QPoint mousePosition();
    static QPoint mouseDelta();
    static void setWheelDelta(QPoint delta);
    static QPoint wheelDelta();

    // State updating
    static void update();
    static void registerKeyPress(int key);
    static void registerKeyRelease(int key);
    static void registerMousePress(Qt::MouseButton button);
    static void registerMouseRelease(Qt::MouseButton button);
    static void reset();
};

inline bool InputManager::keyTriggered(Qt::Key key){
    return keyState(key) == InputTriggered;
}

inline bool InputManager::keyPressed(Qt::Key key){
    return keyState(key) == InputPressed;
}

inline bool InputManager::keyReleased(Qt::Key key){
    return keyState(key) == InputReleased;
}

inline bool InputManager::buttonTriggered(Qt::MouseButton button){
    return buttonState(button) == InputTriggered;
}

inline bool InputManager::buttonPressed(Qt::MouseButton button){
    return buttonState(button) == InputPressed;
}

inline bool InputManager::buttonReleased(Qt::MouseButton button){
    return buttonState(button) == InputReleased;
}

#endif // INPUTMANAGER_H
