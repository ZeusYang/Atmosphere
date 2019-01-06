#ifndef OPENGLTEXTURE_H
#define OPENGLTEXTURE_H

#include "Utilities/OpenGLVersion.h"

class OpenGLTexture{
public:
    OpenGLTexture(GLuint id);
    OpenGLTexture(int width,int height);
    OpenGLTexture(int width,int height,int depth,GLenum format,bool half_precision);
    ~OpenGLTexture();

    void bind(int texture_unit);
    void bind3D(int texture_unit);
    void unBind();

    GLuint getTextureId()const;

private:
    GLuint textureId;
};

#define F(c) c *f = (GlobalContext::contextFunc)
inline void OpenGLTexture::bind(int texture_unit){
    F(OGL_Function);
    f->glActiveTexture(GL_TEXTURE0 + texture_unit);
    f->glBindTexture(GL_TEXTURE_2D, textureId);
}

inline void OpenGLTexture::bind3D(int texture_unit){
    F(OGL_Function);
    f->glActiveTexture(GL_TEXTURE0 + texture_unit);
    f->glBindTexture(GL_TEXTURE_3D, textureId);
}

inline void OpenGLTexture::unBind(){
    F(OGL_Function);
}

inline GLuint OpenGLTexture::getTextureId() const{
    return textureId;
}

#endif // OPENGLTEXTURE_H
