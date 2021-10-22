#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include <stack>
#include "Simulator.h"

#include "tests.h"

#define SIMU_INTERPOLT

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

extern const QColor DefaultColors[];
const uint32_t pathChartCount=DIM_COUNT*(DIM_COUNT-1)/2;
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
                      const BodyVector &,
                      bool * noCollide=nullptr);

    void drawConservativeCharts();
    void drawPathCharts();

private slots:
    void on_timeSlider_valueChanged(int value);

    void on_timeSlider_sliderPressed();

    void on_timeSlider_sliderReleased();

    void on_cbShowCurves_stateChanged(int arg1);

    void on_cbShowDots_stateChanged(int arg1);

private:
    Ui::MainWindow *ui;
    Simulator Simu;
    QChartView * MotionView;
    QChartView * EnergyView;
    std::array<QChartView * ,pathChartCount> pathViews;
    std::array<std::pair<uint16_t,uint16_t>,pathChartCount> dimPairs;

    static QChart * createEmptyChart();

    static void addSeriesToChart(QChart *,QList<QVector<QPointF>>&,
                                 int * =nullptr);
    static void createScatters(QChart*);
    static QPointF moveScatterIndex(QChart*,int index);
};
#endif // MAINWINDOW_H
