#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QWidget>

namespace Ui {
class MainWindow;
}

class DrawCanvas;
class MainWindow : public QMainWindow{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void receiveFPS(QString fps);

    void on_spectrum_ComboBox_currentIndexChanged(const QString &arg1);

    void on_ozone_ComboBox_currentIndexChanged(const QString &arg1);

    void on_float_ComboBox_currentIndexChanged(const QString &arg1);

    void on_texture_ComboBox_currentIndexChanged(const QString &arg1);

    void on_exposure_SpinBox_valueChanged(int arg1);

    void on_num_scattering_SpinBox_valueChanged(int arg1);

    void on_lightshaft_CheckBox_stateChanged(int arg1);

    void on_rayleighCheckBox_stateChanged(int arg1);

    void on_mieCheckBox_stateChanged(int arg1);

    void on_mie_g_DoubleSpinBox_valueChanged(double arg1);

private:
    Ui::MainWindow *ui;
    DrawCanvas *canvas;
};

#endif // MAINWINDOW_H
