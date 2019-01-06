#ifndef DRAWCANVAS_H
#define DRAWCANVAS_H

#include "Utilities/OpenGLWidget.h"
#include "Utilities/Geometry/Quad.h"
#include "Utilities/ShaderManager.h"
#include "Utilities/Entity.h"
#include "Demo/Demo.h"

class ThirdPersonCamera;
class Player;
class DrawCanvas : public OpenGLWidget{
    Q_OBJECT
public:
    DrawCanvas(UpdateBehavior updateBehavior = NoPartialUpdate, QWidget *parent = 0);
    ~DrawCanvas();

    // Setter
    void setExposure(double expos);
    void setMieGCof(double g);
    void setNumscattering(unsigned int nums);
    void setParameters(QString target,bool use);

protected:
    // OpenGL Methods
    void initializeGL()override;
    void resizeGL(int width, int height)override;
    void paintGL()override;
    void updateEvent(DUpdateEvent *event)override;
signals:
    void sendFPS(QString fps);
private:
    Demo *demo;
    bool m_dirty;
    bool tag_spectrum_use;
    bool tag_ozone_use;
    bool tag_halfpre_use;
    bool tag_texture_use;
};

#endif // DRAWCANVAS_H
