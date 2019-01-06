#include "Sphere.h"

Sphere::Sphere(int m, int n){
    longitude   = m;
    latitude    = n;
    vertexCount = (m - 1) * n + 2;
    indexCount  = (m - 2) * n * 6 + 2 * n * 3;

    vertices.resize(vertexCount);
    normals.resize(vertexCount);
    texcoords.resize(vertexCount);
    indices.resize(indexCount);
    initFaces();
}

Sphere::Sphere(const Sphere &rhs){
    vertexCount = rhs.vertexCount;
    indexCount  = rhs.indexCount;

    vertices.assign(rhs.vertices.begin(),rhs.vertices.end());
    normals.assign(rhs.normals.begin(),rhs.normals.end());
    texcoords.assign(rhs.texcoords.begin(),rhs.texcoords.end());
    indices.assign(rhs.indices.begin(),rhs.indices.end());
}

void Sphere::initFaces(){
    float stepAngZ  = M_PI / longitude;
    float stepAngXY = M_PI * 2.0f / latitude;
    uint current = 0;

    float x0 = 0; float y0 = 0; float z0 = 1;
    float v0 = 0; float u0 = 0;
    vertices[current]   = QVector4D(x0,y0,z0,1.0);
    normals[current]    = QVector3D(x0,y0,z0);
    texcoords[current]  = QVector2D(u0,v0);
    ++current;

    // Calculate the vertices according to the equation of sphere
    for (int i = 1; i < longitude; i++) {
        for (int j = 0; j < latitude; j++) {
            float x           = sin(i * stepAngZ) * cos(j * stepAngXY);
            float y           = sin(i * stepAngZ) * sin(j * stepAngXY);
            float z           = cos(i * stepAngZ);
            float v           = (i * stepAngZ) / M_PI;
            float u           = (j * stepAngXY) / M_PI*2.0f;
            vertices[current] = QVector4D(x,y,z,1.0);
            normals[current]  = QVector3D(x,y,z);
            texcoords[current]= QVector2D(u,v);
            ++current;
        }
    }

    z0 = -1;
    v0 = 1; u0 = 0;
    vertices[current]   = QVector4D(x0,y0,z0,1.0);
    normals[current]    = QVector3D(x0,y0,z0);
    texcoords[current]  = QVector2D(u0,v0);

    // Calculate the indices of sphere
    uint curIndex = 0, baseVertex = 0;
    // The first lap
    for (int i = 0; i < latitude; i++) {
        indices[curIndex++] = 0;
        indices[curIndex++] = i + 1;
        if (i != latitude - 1)indices[curIndex] = i + 2;
        else indices[curIndex] = 1;
        ++curIndex;
    }
    // The middle laps
    baseVertex = 1;
    for (int i = 1; i < longitude - 1; i++) {
        for (int j = 0; j < latitude; j++) {
            indices[curIndex++] = baseVertex + j;
            indices[curIndex++] = baseVertex + latitude + j;
            if (j != latitude - 1)
                indices[curIndex] = baseVertex + latitude + (j + 1);
            else
                indices[curIndex] = baseVertex + latitude;
            ++curIndex;

            if (j != latitude - 1) {
                indices[curIndex++] = baseVertex + latitude + (j + 1);
                indices[curIndex++] = baseVertex + (j + 1);
            } else {
                indices[curIndex++] = baseVertex + latitude;
                indices[curIndex++] = baseVertex;
            }
            indices[curIndex++] = baseVertex + j;
        }
        baseVertex += latitude;
    }
    // The last lap
    for (int i = 0; i < latitude; i++) {
        indices[curIndex++] = baseVertex + i;
        indices[curIndex++] = baseVertex + latitude;
        if (i != latitude - 1)indices[curIndex] = baseVertex + i + 1;
        else indices[curIndex] = baseVertex;
        ++curIndex;
    }
}

Sphere::~Sphere() {}
