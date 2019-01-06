#ifndef OPENGLVERSION_H
#define OPENGLVERSION_H
#include <QOpenGLFunctions_4_3_Core>

/***************************
 * OpenGL版本
 ***************************/
namespace GlobalContext{

extern QOpenGLFunctions_4_3_Core* contextFunc;

}

typedef QOpenGLFunctions_4_3_Core OGL_Function;

#endif // OPENGLVERSION_H
