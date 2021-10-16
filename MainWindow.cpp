#include "MainWindow.h"
#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow) {
    ui->setupUi(this);

    ui->conervativeLayout->addWidget(
                conservativeQuatity[0].first=new QChartView(this)
            ,0,0);
    conservativeQuatity[0].second=new QChart();
    conservativeQuatity[0].second->setTitle("Motion");
    conservativeQuatity[0].second->setTheme(QChart::ChartTheme::ChartThemeLight);
    {
        QLineSeries * temp=new QLineSeries(this);
        temp->clear();
        for(int i=0;i<100;i++) {
            temp->append(i,sin(i));
        }
        conservativeQuatity[0].second->addSeries(temp);

        QValueAxis * tempAxis=new QValueAxis;
        tempAxis->setGridLineVisible(true);
        tempAxis->setVisible(true);
        tempAxis->setRange(0,99);
        tempAxis->setTitleText("Time / s");
        conservativeQuatity[0].second->addAxis(tempAxis,Qt::AlignBottom);

        tempAxis=new QValueAxis;
        tempAxis->setGridLineVisible(true);
        tempAxis->setVisible(true);
        tempAxis->setRange(-1,1);
        tempAxis->setTitleText("Motion(SI)");
        conservativeQuatity[0].second->addAxis(tempAxis,Qt::AlignLeft);
        conservativeQuatity[0].second->legend()->hide();
    }
    conservativeQuatity[0].first->setChart(conservativeQuatity[0].second);
    conservativeQuatity[0].first->setRenderHint(QPainter::Antialiasing);

    ui->conervativeLayout->addWidget(
                conservativeQuatity[1].first=new QChartView(this)
            ,1,0);
    conservativeQuatity[1].second=new QChart();
    conservativeQuatity[1].second->setTitle("Energy");
    conservativeQuatity[1].second->setTheme(QChart::ChartTheme::ChartThemeLight);
    conservativeQuatity[1].first->setChart(conservativeQuatity[1].second);
    conservativeQuatity[1].first->setRenderHint(QPainter::Antialiasing);

}

MainWindow::~MainWindow() {
    delete ui;
}

