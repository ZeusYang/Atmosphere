#ifndef DUPDATEEVENT_H
#define DUPDATEEVENT_H

#include <QEvent>

class DUpdateEvent : public QEvent
{
public:
    DUpdateEvent();
    // Static Access
    static QEvent::Type type();
};

inline DUpdateEvent::DUpdateEvent() : QEvent(DUpdateEvent::type()) {}
#endif // DUPDATEEVENT_H
