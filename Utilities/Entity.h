#ifndef ENTITY_H
#define ENTITY_H
#include "Utilities/Geometry/Mesh.h"
#include "Utilities/DataBuffer.h"
#include <QMatrix4x4>

enum VBO_TYPE{
    VERTEX_VBO          = 0,
    NORMAL_VBO          = 1,
    TEXCOORD_VBO        = 2,
    COLOR_VBO           = 3,
};

enum VBO_LOCATION{
    VERTEX_LOCATION     = 0,
    NORMAL_LOCATION     = 1,
    TEXCOORD_LOCATION   = 2,
    COLOR_LOCATION      = 3,
};
class ShaderProgram;
class Entity{
public:
    Entity(Mesh &mesh);
    ~Entity();

    void draw(ShaderProgram *shader) const;
    void setModelMatrix(QMatrix4x4 &target);

private:
    int vertexCount,indexCount;
    bool indexed;

    QMatrix4x4 modelMatrix;
    QMatrix4x4 normalMatrix;

    DataBuffer *buffer;

    std::vector<float> vertexBuffer;
    std::vector<float> normalBuffer;
    std::vector<float> texcoordBuffer;
};

#endif // ENTITY_H
