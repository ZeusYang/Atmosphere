#include "Entity.h"
#include "Utilities/OpenGLVersion.h"
#include "Utilities/ShaderProgram.h"
#define F(c) c *f = (GlobalContext::contextFunc)

Entity::Entity(Mesh &mesh){
    // Push the data into buffers
    vertexCount = mesh.vertexCount;
    indexCount  = mesh.indexCount;

    vertexBuffer.resize(vertexCount * 3);
    normalBuffer.resize(vertexCount * 3);
    texcoordBuffer.resize(vertexCount * 2);

    int storeVertexCount = 0;
    for(auto x = 0;x < vertexCount;++x){
        vertexBuffer[storeVertexCount*3]    = mesh.vertices[x].x();
        vertexBuffer[storeVertexCount*3+1]  = mesh.vertices[x].y();
        vertexBuffer[storeVertexCount*3+2]  = mesh.vertices[x].z();

        normalBuffer[storeVertexCount*3]    = mesh.normals[x].x();
        normalBuffer[storeVertexCount*3+1]  = mesh.normals[x].y();
        normalBuffer[storeVertexCount*3+2]  = mesh.normals[x].z();

        texcoordBuffer[storeVertexCount*2]  = mesh.texcoords[x].x();
        texcoordBuffer[storeVertexCount*2+1]= mesh.texcoords[x].y();

        ++storeVertexCount;
    }
    indexed = !mesh.indices.empty();
    int bufCount = 3;
    bufCount = indexed?bufCount+1:bufCount;

    buffer = new DataBuffer(bufCount);
    buffer->pushData(0, new Data(VERTEX_LOCATION, GL_FLOAT, vertexCount, 3, 1,
                                 buffer->vbos[0], false, GL_STATIC_DRAW, -1, &vertexBuffer[0]));
    buffer->pushData(1, new Data(NORMAL_LOCATION, GL_FLOAT, vertexCount, 3, 1,
                                 buffer->vbos[1], false, GL_STATIC_DRAW, -1, &normalBuffer[0]));
    buffer->pushData(2, new Data(TEXCOORD_LOCATION, GL_FLOAT, vertexCount, 2, 1,
                                 buffer->vbos[2], false, GL_STATIC_DRAW, -1, &texcoordBuffer[0]));

    // Index
    if(indexed){
        buffer->pushData(3, new Data(GL_UNSIGNED_INT, indexCount, buffer->vbos[3],
                         GL_STATIC_DRAW, &mesh.indices[0]));
    }

    buffer->release();
    modelMatrix.setToIdentity();
    normalMatrix.setToIdentity();
}

Entity::~Entity(){
    if(buffer)delete buffer;
    buffer = nullptr;
}

void Entity::draw(ShaderProgram *shader) const{
    // Draw the entity
    F(OGL_Function);
    Q_UNUSED(shader);
    //shader->setMatrix4("modelMatrix",modelMatrix);
    buffer->use();
    if(!indexed){
        f->glDrawArrays(GL_TRIANGLES,0,vertexCount);
    }
    else{
        f->glDrawElements(GL_TRIANGLES,indexCount,GL_UNSIGNED_INT,0);
    }
    buffer->release();
}

void Entity::setModelMatrix(QMatrix4x4 &target){
    this->modelMatrix = target;
}
