#include "OpenGLWidget.h"
#include "Utilities/InputManager.h"
#include "Utilities/OpenGLVersion.h"
#include "Utilities/DUpdateEvent.h"
#include "Utilities/OpenGLFrameTimer.h"

#include <QTimer>
#include <QKeyEvent>
#include <QApplication>
#include <QOpenGLDebugLogger>
#include <QOpenGLDebugMessage>

/**********************************
 * OpenGLWindowPrivate
 **********************************/
class OpenGLWidgetPrivate {
public:
    OpenGLWidgetPrivate(QObject *parent = nullptr);

    // Rendering Statistics
    OpenGLFrameTimer m_frameTimer;
    QString m_fps;
};

OpenGLWidgetPrivate::OpenGLWidgetPrivate(QObject *parent) :
    m_frameTimer(parent),m_fps(""){}


/**********************************
 * OpenGLWidget
 **********************************/

OpenGLWidget::OpenGLWidget(UpdateBehavior updateBehavior, QWidget *parent):
    QOpenGLWidget(parent),m_private(new OpenGLWidgetPrivate(this)){
    this->setUpdateBehavior(updateBehavior);
    // Display fps
    connect(&m_private->m_frameTimer, SIGNAL(timeout(float)),
            this, SLOT(frameTimeout(float)));
}

OpenGLWidget::~OpenGLWidget(){
    makeCurrent();
    delete m_private;
}

/**********************************
 * OpenGLWidget Protected Method
 **********************************/
void OpenGLWidget::initializeGL(){

    // Initialize
    initializeOpenGLFunctions();
    GlobalContext::contextFunc = dynamic_cast<QOpenGLFunctions_4_3_Core*>(this);
    connect(this, SIGNAL(frameSwapped()), this, SLOT(update()));
    connect(this, SIGNAL(frameSwapped()), &m_private->m_frameTimer, SLOT(frameSwapped()));

    QOpenGLWidget::initializeGL();
}

void OpenGLWidget::resizeGL(int width, int height){
    QOpenGLWidget::resizeGL(width, height);
}

void OpenGLWidget::paintGL(){
    QOpenGLWidget::paintGL();
}

QString OpenGLWidget::getFPS() const{
    return m_private->m_fps;
}

/**********************************
 * Event Methods
 **********************************/
bool OpenGLWidget::event(QEvent *event){
    if (event->type() == DUpdateEvent::type()){
        updateEvent(static_cast<DUpdateEvent*>(event));
        return true;
    }
    return QOpenGLWidget::event(event);
}

void OpenGLWidget::keyPressEvent(QKeyEvent *event){
    if (event->isAutoRepeat()){
        event->ignore();
    }
    else{
        InputManager::registerKeyPress(event->key());
    }
}

void OpenGLWidget::keyReleaseEvent(QKeyEvent *event){
    if (event->isAutoRepeat()){
        event->ignore();
    }
    else{
        InputManager::registerKeyRelease(event->key());
    }
}

void OpenGLWidget::mousePressEvent(QMouseEvent *event){
    InputManager::registerMousePress(event->button());
}

void OpenGLWidget::mouseReleaseEvent(QMouseEvent *event){
    InputManager::registerMouseRelease(event->button());
}

void OpenGLWidget::wheelEvent(QWheelEvent *event){
    InputManager::setWheelDelta(event->angleDelta()/15);
}

void OpenGLWidget::moveEvent(QMoveEvent *event){
    QOpenGLWidget::moveEvent(event);
}

void OpenGLWidget::updateEvent(DUpdateEvent *event){
    Q_UNUSED(event);
}

/**********************************
 * Public Slots
 **********************************/
void OpenGLWidget::update(){
    InputManager::update();
    // Update Logic
    {
        DUpdateEvent ev;
        QApplication::sendEvent(this, &ev);
    }

    QOpenGLWidget::update();
}

void OpenGLWidget::frameTimeout(float fps){
    // Frame rate per second
    m_private->m_fps = tr("%1").arg(std::round(fps));
}

/**********************************
 * Public Methods
 **********************************/
void OpenGLWidget::printVersionInformation(){
    QString glType;
    QString glVersion;
    QString glProfile;

    // Get Version Information
    glType = (context()->isOpenGLES()) ? "OpenGL ES" : "OpenGL";
    glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));

    // Get Profile Information
#define CASE(c) case QSurfaceFormat::c: glProfile = #c; break
    switch (format().profile())
    {
    CASE(NoProfile);
    CASE(CoreProfile);
    CASE(CompatibilityProfile);
    }
#undef CASE

    qDebug() << qPrintable(glType)
             << qPrintable(glVersion)
             << "(" << qPrintable(glProfile) << ")";
}


