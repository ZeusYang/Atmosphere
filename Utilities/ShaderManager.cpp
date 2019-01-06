#include "ShaderManager.h"
#include <QDebug>

ShaderManager::ShaderManager(){}

ShaderManager::~ShaderManager(){
    for(auto it = shaders.begin();it != shaders.end();++it)
        if(it->second) delete it->second;
    shaders.clear();
}

ShaderProgram *ShaderManager::addShader(const QString name, const QString comp){
    ShaderProgram* shader = new ShaderProgram(comp);
    shader->name = name;
    shaders.insert(std::pair<QString,ShaderProgram*>(name,shader));
    return shader;
}

ShaderProgram *ShaderManager::addShader(const QString name, const QString vert,
                                        const QString frag, const QString geom){
    ShaderProgram* shader = new ShaderProgram(vert,frag,geom);
    shader->name = name;
    shaders.insert(std::pair<QString,ShaderProgram*>(name,shader));
    return shader;
}

ShaderProgram *ShaderManager::findShader(const QString name)
{
    auto it = shaders.find(name);
    if(it != shaders.end())return it->second;
    return nullptr;
}
