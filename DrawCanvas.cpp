#include "DrawCanvas.h"
#include <QDebug>
#include "Utilities/Player.h"
#include "Utilities/Geometry/Cube.h"
#include "Utilities/Geometry/Sphere.h"
#include "Utilities/ThirdPersonCamera.h"
#include "Utilities/OpenGLVersion.h"

#define F(c) c *f = (GlobalContext::contextFunc)

DrawCanvas::DrawCanvas(UpdateBehavior updateBehavior, QWidget *parent) :
    OpenGLWidget(updateBehavior, parent) ,demo(nullptr),tag_spectrum_use(false),
    tag_ozone_use(true),tag_halfpre_use(true),tag_texture_use(true),m_dirty(false){}

DrawCanvas::~DrawCanvas(){
    makeCurrent();
    if(demo)delete demo;
}

void DrawCanvas::setExposure(double expos){
    demo->setExposure(expos);
}

void DrawCanvas::setMieGCof(double g){
    m_dirty = demo->setMieGCof(g);
}

void DrawCanvas::setNumscattering(unsigned int nums){
    m_dirty = demo->setNumScattering(static_cast<unsigned int>(nums));
}

void DrawCanvas::setParameters(QString target, bool use){
    if(target == QString("Solar spectrum")){
        m_dirty = demo->setSpectrum(use);
    }
    else if(target == QString("Ozone")){
        m_dirty = demo->setOZone(use);
    }
    else if(target == QString("Half precision")){
        m_dirty = demo->setHalfpre(use);
    }
    else if(target == QString("Combined texture")){
        m_dirty = demo->setCombined(use);
    }
    else if(target == QString("Light shaft")){
        demo->setLightshaft(use);
    }
    else if(target == QString("Rayleigh")){
        m_dirty = demo->setRayleigh(use);
    }
    else if(target == QString("Mie")){
        m_dirty = demo->setMie(use);
    }
}

/**********************************
 * OpenGL Methods
 **********************************/
void DrawCanvas::initializeGL(){
    OpenGLWidget::initializeGL();
    printVersionInformation();
    demo = new Demo();
}

void DrawCanvas::resizeGL(int width, int height){
    OpenGLWidget::resizeGL(width, height);
    F(OGL_Function);
    f->glViewport(0,0,this->width(),this->height());
    demo->resize(width,height);
}

void DrawCanvas::paintGL(){
    // 更新模型
    if(m_dirty){
        demo->initModel();
        m_dirty = false;
    }
    F(OGL_Function);

    // Set rendering state
    f->glViewport(0,0,width(),height());
    f->glClearColor(0.5,0.5,0.5,1.0);
    f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    f->glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    f->glEnable(GL_DEPTH_TEST);
    f->glDisable(GL_CULL_FACE);

    // Demo rendering
    demo->draw();

    OpenGLWidget::paintGL();
}

/**********************************
 * Events
 **********************************/
void DrawCanvas::updateEvent(DUpdateEvent *event){
    Q_UNUSED(event);
    demo->update();
    emit sendFPS(getFPS());
}
