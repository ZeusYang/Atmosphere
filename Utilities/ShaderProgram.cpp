#include "ShaderProgram.h"
#include <QDebug>
#include <QFile>
#include <iostream>
#include "Utilities/OpenGLVersion.h"
#include <QOpenGLShaderProgram>

#define F(c) c *f = (GlobalContext::contextFunc)

void ShaderProgram::linkShader(std::vector<uint> shaders){
    F(OGL_Function);
    if(!program->create()){
        qDebug() << "Program Creation Error!";
    }
    for(auto it = shaders.begin();it != shaders.end();++it){
        qDebug() << *it;
        f->glAttachShader(program->programId(),*it);
    }
    f->glLinkProgram(program->programId());
    for(auto it = shaders.begin();it != shaders.end();++it){
        f->glDetachShader(program->programId(),*it);
    }
}

ShaderProgram::ShaderProgram(){
    program = new QOpenGLShaderProgram();
}

ShaderProgram::ShaderProgram(const QString &comp){
    program = new QOpenGLShaderProgram;
    program->addShaderFromSourceFile(QOpenGLShader::Compute,comp);
    program->link();
}

ShaderProgram::ShaderProgram(const QString &vert, const QString &frag,
                             const QString &geom){
    program = new QOpenGLShaderProgram;
    program->addShaderFromSourceFile(QOpenGLShader::Vertex,vert);
    program->addShaderFromSourceFile(QOpenGLShader::Fragment,frag);
    if(!geom.isEmpty()){
        program->addShaderFromSourceFile(QOpenGLShader::Geometry,geom);
    }
    QString logger = program->log();
    if(!logger.isEmpty()){
        qDebug() << logger;
    }
    program->link();
    QString logger_lnk = program->log();
    if(!logger_lnk.isEmpty()){
        qDebug() << logger_lnk;
    }
}

ShaderProgram::ShaderProgram(const char *vert, const char *frag, const char *geom){
    program = new QOpenGLShaderProgram;
    program->addShaderFromSourceCode(QOpenGLShader::Vertex,vert);
    program->addShaderFromSourceCode(QOpenGLShader::Fragment,frag);
    if(geom){
        program->addShaderFromSourceCode(QOpenGLShader::Geometry,geom);
    }
    QString logger = program->log();
    if(!logger.isEmpty()){
        qDebug() << logger;
    }
    program->link();
    QString logger_lnk = program->log();
    if(!logger_lnk.isEmpty()){
        qDebug() << logger_lnk;
    }
}

ShaderProgram::~ShaderProgram(){
    if(program)delete program;
    program = nullptr;
    uniformLocation.clear();
}

void ShaderProgram::use(){
    if(this == nullptr || program == nullptr){
        qDebug() << "a nullptr shader program";
    }
    program->bind();
}

void ShaderProgram::release(){
    program->release();
}

void ShaderProgram::addUniform(const QString &target){
    // Add uniform location that haven't been added before.
    GLuint location = program->uniformLocation(target);
    if(static_cast<int>(location) != INVALID_LOCATION){
        uniformLocation.insert(std::pair<QString,uint>(target,location));
    }
}

void ShaderProgram::setInt(const QString &param, int value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
    else qDebug() << "INVALID_LOCATION->"  << param;
}

void ShaderProgram::setFloat(const QString &param, float value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
    else qDebug() << "INVALID_LOCATION->"  << param;
}

void ShaderProgram::setVector2(const QString &param, float x, float y){
    setVector2(param,QVector2D(x,y));
}

void ShaderProgram::setVector2(const QString &param, QVector2D value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
    else qDebug() << "INVALID_LOCATION->"  << param;
}

void ShaderProgram::setVector3(const QString &param, float x, float y, float z){
    setVector3(param,QVector3D(x,y,z));
}

void ShaderProgram::setVector3(const QString &param, QVector3D value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
    else qDebug() << "INVALID_LOCATION->"  << param;
}

void ShaderProgram::setVector4(const QString &param,float x,float y,float z,float w){
    setVector4(param,QVector4D(x,y,z,w));
}

void ShaderProgram::setVector4(const QString &param, QVector4D value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
    else qDebug() << "INVALID_LOCATION->"  << param;
}

void ShaderProgram::setMatrix2(const QString &param, QMatrix2x2 value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
}

void ShaderProgram::setMatrix3(const QString &param, QMatrix3x3 value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
}

void ShaderProgram::setMatrix4(const QString &param, QMatrix4x4 value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if (location != INVALID_LOCATION)program->setUniformValue(location,value);
    else qDebug() << "INVALID_LOCATION->"  << param;
}

void ShaderProgram::setMatrix4(const QString &param, int count,
                               std::vector<QMatrix4x4> &value){
    int location = findUniformLocations(param);
    if (location == INVALID_LOCATION) {
        addUniform(param);
        location = findUniformLocations(param);
    }
    if(location != INVALID_LOCATION)
        program->setUniformValueArray(location,&value[0],count);
    else qDebug() << "INVALID_LOCATION->" << param;
}

QString ShaderProgram::getShaderFromFile(QString path) {
    QFile file(path);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QByteArray con = file.readAll();
    QString context(con);
    file.close();
    return context;
}

int ShaderProgram::findUniformLocations(const QString &target){
    auto it = uniformLocation.find(target);
    if(it != uniformLocation.end())
        return it->second;
    return INVALID_LOCATION;
}

