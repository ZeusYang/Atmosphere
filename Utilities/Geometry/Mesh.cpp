#include "Mesh.h"

Mesh::Mesh():
    vertexCount(0),indexCount(0) {}

Mesh::Mesh(const Mesh &rhs){
    Q_UNUSED(rhs);
}

void Mesh::ReleaseData()
{
    vertices.clear();
    normals.clear();
    texcoords.clear();
    indices.clear();
}

Mesh::~Mesh(){}
