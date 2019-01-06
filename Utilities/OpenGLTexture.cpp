#include "OpenGLTexture.h"

#define F(c) c *f = (GlobalContext::contextFunc)

OpenGLTexture::OpenGLTexture(GLuint id)
    :textureId(id) {}

OpenGLTexture::OpenGLTexture(int width, int height){
    F(OGL_Function);
    f->glGenTextures(1, &textureId);
    f->glActiveTexture(GL_TEXTURE0);
    f->glBindTexture(GL_TEXTURE_2D, textureId);
    f->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    f->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    f->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    f->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    f->glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    // 16F precision for the transmittance gives artifacts.
    f->glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0,
                    GL_RGBA, GL_FLOAT, NULL);
}

OpenGLTexture::OpenGLTexture(int width, int height, int depth,
                             GLenum format, bool half_precision){
    F(OGL_Function);
    f->glGenTextures(1, &textureId);
    f->glActiveTexture(GL_TEXTURE0);
    f->glBindTexture(GL_TEXTURE_3D, textureId);
    f->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    f->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    f->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    f->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    f->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    f->glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    GLenum internal_format = format == GL_RGBA ?
                (half_precision ? GL_RGBA16F : GL_RGBA32F) :
                (half_precision ? GL_RGB16F : GL_RGB32F);
    f->glTexImage3D(GL_TEXTURE_3D, 0, internal_format, width, height, depth, 0,
                    format, GL_FLOAT, NULL);
}

OpenGLTexture::~OpenGLTexture(){
    F(OGL_Function);
    f->glDeleteTextures(1, &textureId);
}
