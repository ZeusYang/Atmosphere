#ifndef SHADERPROGRAM_H
#define SHADERPROGRAM_H
#include <map>
#include <QString>
#include <QVector2D>
#include <QVector3D>
#include <QVector4D>
#include <QMatrix2x2>
#include <QMatrix3x3>
#include <QMatrix4x4>

#ifndef INVALID_LOCATION
#define INVALID_LOCATION -1
#endif

class QOpenGLShaderProgram;
class QOpenGLShader;
class ShaderProgram {
public:

    QString name;

    QOpenGLShaderProgram *program;

    void linkShader(std::vector<uint> shaders);

    ShaderProgram();

    ShaderProgram(const QString &comp);

    ShaderProgram(const QString &vert,const QString &frag,const QString &geom = "");

    ShaderProgram(const char *vert,const char *frag,const char *geom = nullptr);

    ~ShaderProgram();

    void use();

    void release();

    void addUniform(const QString &target);

    void     setInt(const QString &param,int value);
    void   setFloat(const QString &param,float value);
    void setVector2(const QString &param,float x,float y);
    void setVector2(const QString &param,QVector2D value);
    void setVector3(const QString &param,float x,float y,float z);
    void setVector3(const QString &param,QVector3D value);
    void setVector4(const QString &param,float x,float y,float z,float w);
    void setVector4(const QString &param,QVector4D value);
    void setMatrix2(const QString &param,QMatrix2x2 value);
    void setMatrix3(const QString &param,QMatrix3x3 value);
    void setMatrix4(const QString &param,QMatrix4x4 value);
    void setMatrix4(const QString &param,int count,std::vector<QMatrix4x4> &value);

    static QString getShaderFromFile(QString path);

private:
    std::map<QString, uint> uniformLocation;

    int findUniformLocations(const QString &target);
};

#endif // SHADERPROGRAM_H
