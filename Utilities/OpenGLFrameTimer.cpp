#include "OpenGLFrameTimer.h"
#include <QElapsedTimer>

/**********************************
 * OpenGLFrameTimerPrivate
 **********************************/
class OpenGLFrameTimerPrivate {
public:
    OpenGLFrameTimerPrivate();

    // Public Members
    int m_frameCount;
    int m_frameCountFrequency;
    QElapsedTimer m_frameTimer;
};

OpenGLFrameTimerPrivate::OpenGLFrameTimerPrivate() :
    m_frameCount(0), m_frameCountFrequency(60){
    m_frameTimer.start();
}

/**********************************
 * OpenGLFrameTimer
 **********************************/
#define P(c) c &p = *m_private

OpenGLFrameTimer::OpenGLFrameTimer(QObject *parent) : QObject(parent)
  ,m_private(new OpenGLFrameTimerPrivate){
}

OpenGLFrameTimer::~OpenGLFrameTimer(){
    delete m_private;
}

/**********************************
 * Public Methods
 **********************************/
void OpenGLFrameTimer::setFrequency(int hz){
    if (hz <= 0) hz = 1;
    m_private->m_frameCountFrequency = hz;
}

int OpenGLFrameTimer::frequency() const{
    return m_private->m_frameCountFrequency;
}

/**********************************
 * Public Slots
 **********************************/
void OpenGLFrameTimer::frameSwapped(){
    //并非真正的帧缓冲交换，仅仅是计算fps，发送信号
    //为了避免频繁计算,设定一个阈值m_frameCountFrequency,仅当达到阈值才计算
    if (++m_private->m_frameCount > m_private->m_frameCountFrequency){
        qint64 ms = m_private->m_frameTimer.elapsed();
        float sec = ms / 1000.0f;
        emit timeout(m_private->m_frameCount / sec);
        m_private->m_frameCount = 0;
        m_private->m_frameTimer.start();
    }
}
