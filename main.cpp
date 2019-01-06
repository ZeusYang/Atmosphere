#include "MainWindow.h"
#include <QApplication>

int main(int argc, char *argv[]){
    QApplication app(argc, argv);

    MainWindow window;
    window.setWindowTitle(QObject::tr("基于GPU预计算的大气层光效渲染"));
    window.show();

    return app.exec();
}
