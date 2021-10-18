#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include <stack>
#include "Simulator.h"

#include "tests.h"

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

    void drawConservativeCharts();
    void drawPathCharts();

private:
    Ui::MainWindow *ui;
    Simulator Simu;
    QChartView * MotionView;
    QChartView * EnergyView;
    std::array<QChartView * ,DIM_COUNT*(DIM_COUNT)/2> pathViews;
    std::array<std::pair<uint16_t,uint16_t>,DIM_COUNT*(DIM_COUNT)/2> dimPairs;

    static QChart * createEmptyChart();

    static void addSeriesToChart(QChart *,QList<QVector<QPointF>>&,
                                 int * =nullptr);

};
#endif // MAINWINDOW_H
