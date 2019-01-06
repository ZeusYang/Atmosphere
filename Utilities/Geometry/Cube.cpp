#include "Cube.h"

Cube::Cube() : Mesh(){
    indexCount  = 36;
    vertexCount = 24;

    vertices.resize(vertexCount);
    normals.resize(vertexCount);
    texcoords.resize(vertexCount);
    indices.resize(indexCount);
    initFaces();
}

Cube::Cube(const Cube &rhs) {
    indexCount  = rhs.indexCount;
    vertexCount = rhs.vertexCount;

    vertices.assign(rhs.vertices.begin(),rhs.vertices.end());
    normals.assign(rhs.normals.begin(),rhs.normals.end());
    texcoords.assign(rhs.texcoords.begin(),rhs.texcoords.end());
    indices.assign(rhs.indices.begin(),rhs.indices.end());
}

Cube::~Cube(){}

// Front Verticies
#define VERTEX_FTR QVector4D( 0.5f,  0.5f,  0.5f, 1.0f)
#define VERTEX_FTL QVector4D(-0.5f,  0.5f,  0.5f, 1.0f)
#define VERTEX_FBL QVector4D(-0.5f, -0.5f,  0.5f, 1.0f)
#define VERTEX_FBR QVector4D( 0.5f, -0.5f,  0.5f, 1.0f)

// Back Verticies
#define VERTEX_BTR QVector4D( 0.5f,  0.5f, -0.5f, 1.0f)
#define VERTEX_BTL QVector4D(-0.5f,  0.5f, -0.5f, 1.0f)
#define VERTEX_BBL QVector4D(-0.5f, -0.5f, -0.5f, 1.0f)
#define VERTEX_BBR QVector4D( 0.5f, -0.5f, -0.5f, 1.0f)

void Cube::initFaces(){
    //front face
    vertices[0]     = VERTEX_FTR;
    vertices[1]     = VERTEX_FTL;
    vertices[2]     = VERTEX_FBL;
    vertices[3]     = VERTEX_FBR;
    normals[0]      = QVector3D(0,0,1);
    normals[1]      = QVector3D(0,0,1);
    normals[2]      = QVector3D(0,0,1);
    normals[3]      = QVector3D(0,0,1);
    texcoords[0]    = QVector2D(1,1);
    texcoords[1]    = QVector2D(0,1);
    texcoords[2]    = QVector2D(0,0);
    texcoords[3]    = QVector2D(1,0);
    indices[0]      = 0;indices[1] = 1;indices[2] = 2;
    indices[3]      = 0;indices[4] = 2;indices[5] = 3;

    //back face
    vertices[4]     = VERTEX_BTL;
    vertices[5]     = VERTEX_BTR;
    vertices[6]     = VERTEX_BBR;
    vertices[7]     = VERTEX_BBL;
    normals[4]      = QVector3D(0,0,-1);
    normals[5]      = QVector3D(0,0,-1);
    normals[6]      = QVector3D(0,0,-1);
    normals[7]      = QVector3D(0,0,-1);
    texcoords[4]    = QVector2D(1,1);
    texcoords[5]    = QVector2D(0,1);
    texcoords[6]    = QVector2D(0,0);
    texcoords[7]    = QVector2D(1,0);
    indices[6]      = 4;indices[7] = 5;indices[8] = 6;
    indices[9]      = 4;indices[10] = 6;indices[11] = 7;

    //top face
    vertices[8]     = VERTEX_BTR;
    vertices[9]     = VERTEX_BTL;
    vertices[10]    = VERTEX_FTL;
    vertices[11]    = VERTEX_FTR;
    normals[8]      = QVector3D(0,1,0);
    normals[9]      = QVector3D(0,1,0);
    normals[10]     = QVector3D(0,1,0);
    normals[11]     = QVector3D(0,1,0);
    texcoords[8]    = QVector2D(1,1);
    texcoords[9]    = QVector2D(0,1);
    texcoords[10]   = QVector2D(0,0);
    texcoords[11]   = QVector2D(1,0);
    indices[12]     = 8;indices[13] = 9;indices[14] = 10;
    indices[15]     = 8;indices[16] = 10;indices[17] = 11;

    //bottom face
    vertices[12]    = VERTEX_BBL;
    vertices[13]    = VERTEX_BBR;
    vertices[14]    = VERTEX_FBR;
    vertices[15]    = VERTEX_FBL;
    normals[12]     = QVector3D(0,-1,0);
    normals[13]     = QVector3D(0,-1,0);
    normals[14]     = QVector3D(0,-1,0);
    normals[15]     = QVector3D(0,-1,0);
    texcoords[12]   = QVector2D(1,1);
    texcoords[13]   = QVector2D(0,1);
    texcoords[14]   = QVector2D(0,0);
    texcoords[15]   = QVector2D(1,0);
    indices[18]     = 12;indices[19] = 13;indices[20] = 14;
    indices[21]     = 12;indices[22] = 14;indices[23] = 15;

    //left face
    vertices[16]    = VERTEX_FTL;
    vertices[17]    = VERTEX_BTL;
    vertices[18]    = VERTEX_BBL;
    vertices[19]    = VERTEX_FBL;
    normals[16]     = QVector3D(-1,0,0);
    normals[17]     = QVector3D(-1,0,0);
    normals[18]     = QVector3D(-1,0,0);
    normals[19]     = QVector3D(-1,0,0);
    texcoords[16]   = QVector2D(1,1);
    texcoords[17]   = QVector2D(0,1);
    texcoords[18]   = QVector2D(0,0);
    texcoords[19]   = QVector2D(1,0);
    indices[24]     = 16;indices[25] = 17;indices[26] = 18;
    indices[27]     = 16;indices[28] = 18;indices[29] = 19;

    //right face
    vertices[20]    = VERTEX_BTR;
    vertices[21]    = VERTEX_FTR;
    vertices[22]    = VERTEX_FBR;
    vertices[23]    = VERTEX_BBR;
    normals[20]     = QVector3D(1,0,0);
    normals[21]     = QVector3D(1,0,0);
    normals[22]     = QVector3D(1,0,0);
    normals[23]     = QVector3D(1,0,0);
    texcoords[20]   = QVector2D(1,1);
    texcoords[21]   = QVector2D(0,1);
    texcoords[22]   = QVector2D(0,0);
    texcoords[23]   = QVector2D(1,0);
    indices[30]     = 20;indices[31] = 21;indices[32] = 22;
    indices[33]     = 20;indices[34] = 22;indices[35] = 23;
}

#undef VERTEX_BBR
#undef VERTEX_BBL
#undef VERTEX_BTL
#undef VERTEX_BTR

#undef VERTEX_FBR
#undef VERTEX_FBL
#undef VERTEX_FTL
#undef VERTEX_FTR
