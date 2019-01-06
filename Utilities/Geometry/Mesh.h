#ifndef MESH_H
#define MESH_H

#include <QVector3D>
#include <QVector2D>
#include <QVector4D>
#include <vector>

class Mesh{
public:
    int vertexCount,indexCount;
    std::vector<unsigned int>   indices;        //索引
    std::vector<QVector4D>      vertices;       //顶点
    std::vector<QVector3D>      normals;        //法线
    std::vector<QVector2D>      texcoords;      //UV

    Mesh();
    Mesh(const Mesh& rhs);
    void ReleaseData();
    virtual ~Mesh();

private:
    virtual void initFaces() = 0;
};

#endif // MESH_H
