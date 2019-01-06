#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "DrawCanvas.h"
#include <QHBoxLayout>

static bool checkVersion(QOpenGLContext &context, QSurfaceFormat &format);
static QSurfaceFormat* getFirstSupported(std::vector<QSurfaceFormat> &formats);

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow){
    ui->setupUi(this);
    canvas = new DrawCanvas(QOpenGLWidget::NoPartialUpdate,
                            qobject_cast<QWidget*>(this));
    std::vector<QSurfaceFormat> formats;

    // OpenGL4.3 core version
    QSurfaceFormat glFormat;
    glFormat.setRenderableType(QSurfaceFormat::OpenGL);
    glFormat.setProfile(QSurfaceFormat::CoreProfile);
    glFormat.setVersion(4,3);
    formats.push_back(glFormat);
    // Find out which version we support
    QSurfaceFormat *format = getFirstSupported(formats);
    if (format == NULL){
        qFatal("No valid supported version of OpenGL found on device!");
    }
    format->setDepthBufferSize(8);
    canvas->setFormat(*format);

    // 将canvas附加到widget上
    QHBoxLayout *layout = new QHBoxLayout(ui->panel);
    layout->addWidget(canvas);
    canvas->setFocusPolicy(Qt::StrongFocus);

    // 设置Instructions文本框透明
    QPalette pl = ui->instructions->palette();
    pl.setBrush(QPalette::Base,QBrush(QColor(255,0,0,0)));
    ui->instructions->setPalette(pl);

    // fps
    connect(canvas,&DrawCanvas::sendFPS,this,&MainWindow::receiveFPS);
}

MainWindow::~MainWindow(){
    delete ui;
}

void MainWindow::receiveFPS(QString fps){
    ui->fps_LineEdit->setText(fps);
}

static bool checkVersion(QOpenGLContext &context, QSurfaceFormat &format){
    // Checking whether satifying the requirement
    QSurfaceFormat currSurface  = context.format();
    QPair<int,int> currVersion  = currSurface.version();
    QPair<int,int> reqVersion   = format.version();
    if (currVersion.first > reqVersion.first)return true;
    return (currVersion.first == reqVersion.first
            && currVersion.second >= reqVersion.second);
}

static QSurfaceFormat* getFirstSupported(std::vector<QSurfaceFormat> &formats){
    QOpenGLContext context;
    for (QSurfaceFormat &format : formats){
        context.setFormat(format);
        if (context.create()){
            if (checkVersion(context, format)) return &format;
        }
    }
    return NULL;
}

void MainWindow::on_spectrum_ComboBox_currentIndexChanged(const QString &arg1){
    if(arg1 == QString("真实值"))
        canvas->setParameters("Solar spectrum",false);
    else if(arg1 == QString("常数值"))
        canvas->setParameters("Solar spectrum",true);
}

void MainWindow::on_ozone_ComboBox_currentIndexChanged(const QString &arg1){
    if(arg1 == QString("开启"))
        canvas->setParameters("Ozone",true);
    else if(arg1 == QString("关闭"))
        canvas->setParameters("Ozone",false);
}

void MainWindow::on_float_ComboBox_currentIndexChanged(const QString &arg1){
    if(arg1 == QString("半精度"))
        canvas->setParameters("Half precision",true);
    else if(arg1 == QString("单精度"))
        canvas->setParameters("Half precision",false);
}

void MainWindow::on_texture_ComboBox_currentIndexChanged(const QString &arg1){
    if(arg1 == QString("开启"))
        canvas->setParameters("Combined texture",true);
    else if(arg1 == QString("关闭"))
        canvas->setParameters("Combined texture",false);
}

void MainWindow::on_exposure_SpinBox_valueChanged(int arg1){
    canvas->setExposure(static_cast<double>(arg1));
}

void MainWindow::on_num_scattering_SpinBox_valueChanged(int arg1){
    canvas->setNumscattering(arg1);
}

void MainWindow::on_lightshaft_CheckBox_stateChanged(int arg1){
    if(arg1)canvas->setParameters("Light shaft",true);
    else canvas->setParameters("Light shaft",false);
}

void MainWindow::on_rayleighCheckBox_stateChanged(int arg1){
    if(arg1)canvas->setParameters("Rayleigh",true);
    else canvas->setParameters("Rayleigh",false);
}

void MainWindow::on_mieCheckBox_stateChanged(int arg1){
    if(arg1)canvas->setParameters("Mie",true);
    else canvas->setParameters("Mie",false);
}

void MainWindow::on_mie_g_DoubleSpinBox_valueChanged(double arg1){
    canvas->setMieGCof(arg1);
}
