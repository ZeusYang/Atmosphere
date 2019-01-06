#ifndef SHADERMANAGER_H
#define SHADERMANAGER_H
#include "Utilities/ShaderProgram.h"

class ShaderManager
{
public:
    ShaderManager();

    ~ShaderManager();

    ShaderProgram* addShader(const QString name,const QString comp);

    ShaderProgram* addShader(const QString name,const QString vert,
                             const QString frag,const QString geom = "");

    ShaderProgram* findShader(const QString name);

private:
    std::map<QString,ShaderProgram*> shaders;
};

#endif // SHADERMANAGER_H
