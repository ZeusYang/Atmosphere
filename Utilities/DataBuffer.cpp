#include "DataBuffer.h"
#include "Utilities/OpenGLVersion.h"
/*******************************************
 * Data
 *******************************************/

#define F(c) c *f = (GlobalContext::contextFunc)

Data::Data(uint loc,GLenum type,uint count,uint channel,uint row,const QOpenGLBuffer &vbo,
           bool normalize,GLenum draw,int divisor, void *data)
{
    switch(type){
    case GL_FLOAT:
        bitSize = sizeof(GLfloat);  break;
    case GL_INT:
        bitSize = sizeof(GLint);    break;
    case GL_UNSIGNED_INT:
        bitSize = sizeof(GLuint);   break;
    case GL_UNSIGNED_SHORT:
        bitSize = sizeof(GLushort); break;
    case GL_UNSIGNED_BYTE:
        bitSize = sizeof(GLubyte);  break;
    }
    rowCount     = row;                             //有多少行数据
    channelCount = channel;                         //数据构成数量
    dataSize     = count * channelCount * rowCount; //总数据大小
    bufferid     = vbo;                             //vbo的id
    drawType     = draw;                            //绘制类型,如static_draw
    streamData   = data;                            //数据内存首指针

    //get current ogl context
    F(OGL_Function);
    f->glBindBuffer(GL_ARRAY_BUFFER, bufferid.bufferId());
    f->glBufferData(GL_ARRAY_BUFFER, dataSize * bitSize, streamData, drawType);
    for (uint i = 0; i < row; i++) {
        uint attrloc = loc + i;
        f->glVertexAttribPointer(attrloc, channel, type, normalize,
                                 bitSize * row * channel, (void*)(bitSize * i * channel));
        if (divisor >= 0)
            f->glVertexAttribDivisor(attrloc, divisor);
        f->glEnableVertexAttribArray(attrloc);
    }
    bufferid.release();
}

Data::Data(GLenum type, uint size, const QOpenGLBuffer &vbo, GLenum draw, void *data){
    switch (type) {
    case GL_UNSIGNED_INT:
        bitSize = sizeof(GLuint);break;
    case GL_UNSIGNED_SHORT:
        bitSize = sizeof(GLushort);break;
    case GL_UNSIGNED_BYTE:
        bitSize = sizeof(GLubyte);break;
    }
    dataSize    = size;
    bufferid    = vbo;
    drawType    = draw;
    streamData  = data;
    F(OGL_Function);
    f->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferid.bufferId());
    f->glBufferData(GL_ELEMENT_ARRAY_BUFFER, dataSize * bitSize, streamData, drawType);
}

void Data::updateAttrBuf(uint count, void *data, GLenum draw){
    drawType = draw;
    streamData = data;
    dataSize = count * channelCount * rowCount;
    F(OGL_Function);
    f->glBindBuffer(GL_ARRAY_BUFFER, bufferid.bufferId());
    f->glBufferSubData(GL_ARRAY_BUFFER, 0, dataSize * bitSize, streamData);
}

/*******************************************
 * DataBuffer
 *******************************************/

DataBuffer::DataBuffer(uint nums){
    // Create nums vbos
    bufferNumbers = nums;
    vbos.resize(bufferNumbers);
    // streamDatas stores nums kinds of data
    streamDatas.resize(bufferNumbers);
    for (uint i = 0; i < bufferNumbers; i++) pushData(i, nullptr);
    vao.create();
    vao.bind();
    for(auto &vbo:vbos)vbo.create();
}

void DataBuffer::release(){
    F(OGL_Function);
    f->glBindVertexArray(0);
    f->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    f->glBindBuffer(GL_ARRAY_BUFFER, 0);
}

DataBuffer::~DataBuffer(){
    for(auto &vbo:vbos)vbo.destroy();
    vao.destroy();
    for(auto &elem:streamDatas){
        if(elem)delete elem;
    }
}
