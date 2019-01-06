#ifndef RAWSHADERPROGRAM_H
#define RAWSHADERPROGRAM_H

#include "Utilities/OpenGLVersion.h"
#include <QVector2D>
#include <QVector3D>
#include <QVector4D>
#include <QMatrix4x4>

class RawShaderProgram{
public:
    RawShaderProgram(std::vector<GLuint> shaders);
    ~RawShaderProgram();

    void use()const;
    void release()const;

    GLuint programId();

    void     setInt(const QString &param,int value);
    void   setFloat(const QString &param,float value);
    void setVector2(const QString &param,QVector2D value);
    void setVector3(const QString &param,QVector3D value);
    void setVector4(const QString &param,QVector4D value);
    void setMatrix4(const QString &param,QMatrix4x4 value);
private:
    GLuint id;
};

#endif // RAWSHADERPROGRAM_H
