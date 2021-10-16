#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include <stack>
#include "Simulator.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void runSimulaton(Simulator::Algorithm);

    void runSimulaton(Simulator::Algorithm,
                      const double,
                      const TimeSpan &,
                      const Statue &,
                      const MassVector &,
                      bool * noCollide=nullptr);

    void drawCharts();

private:
    Ui::MainWindow *ui;
    Simulator Simu;
    QChartView * MotionView;
    QChartView * EnergyView;

    static QChart * createEmptyChart();

    static void addSeriesToChart(QChart *,QList<QVector<QPointF>>&,
                                 int * =nullptr);

};
#endif // MAINWINDOW_H
