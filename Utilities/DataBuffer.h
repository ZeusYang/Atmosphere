#ifndef DATABUFFER_H
#define DATABUFFER_H

#include <vector>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>

class Data{
public:
    uint bitSize;
    GLenum drawType;
    void *streamData;
    QOpenGLBuffer bufferid;
    //数据数量、长度
    uint dataSize, channelCount, rowCount;

    Data(uint loc,GLenum type,uint count,uint channel,uint row,
         const QOpenGLBuffer &vbo,bool normalize,GLenum draw,int divisor,
         void *data);
    Data(GLenum type, uint size, const QOpenGLBuffer &vbo, GLenum draw, void* data);
    void updateAttrBuf(uint count, void* data, GLenum draw);
};

class DataBuffer{
public:
    QOpenGLVertexArrayObject vao;
    std::vector<QOpenGLBuffer> vbos;
    std::vector<Data*> streamDatas;
    uint bufferNumbers;             // quantity of vbo

    DataBuffer(uint nums);

    ~DataBuffer();

    void use();                     //like bind

    void release();                 // like unbind

    void pushData(uint index, Data *data);
};

inline void DataBuffer::pushData(uint index, Data *data){
    streamDatas[index] = data;
}

inline void DataBuffer::use(){
    vao.bind();
}

#endif // DATABUFFER_H
