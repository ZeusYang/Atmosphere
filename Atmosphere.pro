#-------------------------------------------------
#
# Project created by QtCreator 2018-12-27T19:05:31
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Atmosphere
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        main.cpp \
        MainWindow.cpp \
    Utilities/OpenGLWidget.cpp \
    Utilities/OpenGLVersion.cpp \
    Utilities/InputManager.cpp \
    Utilities/ShaderProgram.cpp \
    Utilities/ShaderManager.cpp \
    Utilities/DataBuffer.cpp \
    Utilities/Geometry/Mesh.cpp \
    Utilities/Geometry/Quad.cpp \
    Utilities/Geometry/Cube.cpp \
    Utilities/Geometry/Sphere.cpp \
    Utilities/DUpdateEvent.cpp \
    Utilities/OpenGLFrameTimer.cpp \
    DrawCanvas.cpp \
    Utilities/Entity.cpp \
    Utilities/Transform3D.cpp \
    Utilities/Camera3D.cpp \
    Utilities/ThirdPersonCamera.cpp \
    Utilities/Player.cpp \
    Demo/Demo.cpp \
    AtmosphereModel/AtmosphereModel.cpp \
    Utilities/OpenGLTexture.cpp \
    Demo/GodCamera.cpp \
    Demo/RawShaderProgram.cpp

HEADERS += \
        MainWindow.h \
    Utilities/OpenGLWidget.h \
    Utilities/OpenGLVersion.h \
    Utilities/InputManager.h \
    Utilities/ShaderProgram.h \
    Utilities/ShaderManager.h \
    Utilities/DataBuffer.h \
    Utilities/Geometry/Mesh.h \
    Utilities/Geometry/Quad.h \
    Utilities/Geometry/Cube.h \
    Utilities/Geometry/Sphere.h \
    Utilities/DUpdateEvent.h \
    Utilities/OpenGLFrameTimer.h \
    DrawCanvas.h \
    Utilities/Entity.h \
    Utilities/Transform3D.h \
    Utilities/Camera3D.h \
    Utilities/ThirdPersonCamera.h \
    Utilities/Player.h \
    Demo/Demo.h \
    AtmosphereModel/AtmosphereModel.h \
    AtmosphereModel/Constants.h \
    Utilities/OpenGLTexture.h \
    Demo/GodCamera.h \
    Demo/RawShaderProgram.h

FORMS += \
        MainWindow.ui

RESOURCES += \
    resource.qrc
