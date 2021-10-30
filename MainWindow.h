/*
 Copyright © 2021  TokiNoBug
This file is part of ThreeBodySimulation.

    ThreeBodySimulation is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SlopeCraft is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ThreeBodySimulation.  If not, see <https://www.gnu.org/licenses/>.

    Contact with me:
    github:https://github.com/ToKiNoBug
    bilibili:https://space.bilibili.com/351429231
*/

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

    void onUserInput();

    void on_BtnRunSimulation_clicked();

private:
    Ui::MainWindow *ui;
    Simulator Simu;
    QChartView * MotionView;
    QChartView * EnergyView;
    std::array<QChartView * ,pathChartCount> pathViews;
    std::array<std::pair<uint16_t,uint16_t>,pathChartCount> dimPairs;
    std::array<QLineEdit*,BODY_COUNT> massWidgets;
    std::array<std::array<QLineEdit*,BODY_COUNT>,DIM_COUNT> positionWidgets;
    std::array<std::array<QLineEdit*,BODY_COUNT>,DIM_COUNT> velocityWidgets;

    void buildMotionEnergyUI();
    void buildPathUI();
    void buildParameterUI();
    void setParamaters(const BodyVector & mass,
                                    const Statue & y0,
                                    TimeSpan ts,
                                    double step,
                                    Simulator::Algorithm);
    bool grabParameters(BodyVector & mass,
                        Statue & y0,
                        TimeSpan & ts,
                        double & step,
                        Simulator::Algorithm&);

    static QChart * createEmptyChart();

    static void addSeriesToChart(QChart *,QList<QVector<QPointF>>&,
                                 int * =nullptr);
    static void createScatters(QChart*);
    static QPointF moveScatterIndex(QChart*,int index);
};
#endif // MAINWINDOW_H
