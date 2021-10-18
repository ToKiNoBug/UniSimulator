#include "MainWindow.h"
#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow) {
    ui->setupUi(this);


    ui->conervativeLayout->addWidget(
                MotionView=new QChartView(this)
            ,0,0);
    MotionView->setChart(createEmptyChart());
    MotionView->chart()->setTitle("Motion");
    MotionView->chart()->axes(Qt::Horizontal).first()->setTitleText("Time (year)");
    MotionView->chart()->axes(Qt::Vertical).first()->setTitleText("Motion (SI)");
    MotionView->setRenderHint(QPainter::Antialiasing, true);

    ui->conervativeLayout->addWidget(
                EnergyView=new QChartView(this)
            ,1,0);
    EnergyView->setChart(createEmptyChart());
    EnergyView->chart()->setTitle("Energy");
    EnergyView->chart()->axes(Qt::Horizontal).first()->setTitleText("Time (year)");
    EnergyView->chart()->axes(Qt::Vertical).first()->setTitleText("Energy (SI)");
    EnergyView->setRenderHint(QPainter::Antialiasing, true);




    /*
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
*/
}

QChart * MainWindow::createEmptyChart() {
    QChart * chart=new QChart;
    //qDebug()<<__LINE__;
    chart->setTitle("DefaultTitle");
    chart->setTheme(QChart::ChartTheme::ChartThemeLight);
    //qDebug()<<__LINE__;

    QValueAxis * xAxis=new QValueAxis(chart),* yAxis=new QValueAxis(chart);
    xAxis->setRange(0,1);
    yAxis->setRange(0,1);
    xAxis->setTitleText("x");
    yAxis->setTitleText("y");
    xAxis->setGridLineVisible(true);
    yAxis->setGridLineVisible(true);

    //qDebug()<<__LINE__;

    chart->addAxis(xAxis,Qt::AlignBottom);
    chart->addAxis(yAxis,Qt::AlignLeft);

    //qDebug()<<__LINE__;

    chart->layout()->setContentsMargins(5, 5, 5, 5);
    chart->setMargins(QMargins(0, 0, 0, 0));
    QChart::AnimationOptions op;

    op.setFlag(QChart::AnimationOption::GridAxisAnimations,false);
    op.setFlag(QChart::AnimationOption::SeriesAnimations,true);

    chart->setAnimationOptions(op);
    //chart->setBackgroundRoundness(0);

    return chart;
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::runSimulaton(Simulator::Algorithm algo) {
    MassVector mass(Ms,1*Ms);

    Statue start;
    start.first.setValues({{-rs,rs},
                                    {0,0}});
    start.second.setValues({{-0*vs,0},
                                          {-vs,vs}});

    TimeSpan tSpan=std::make_pair(0*year,10*year);

    Time step=0.001*year;

    runSimulaton(algo,step,tSpan,start,mass);
}

void MainWindow::runSimulaton(Simulator::Algorithm algo,
                                                       const double h,
                                                       const TimeSpan & ts,
                                                       const Statue & s,
                                                       const MassVector & mv,
                                                       bool * noCollide) {
    Simulator fast=Simu;

    fast.setMass(mv);

    switch (algo) {
    case Simulator::Algorithm::Euler:
        fast.simulateEuler(h,ts,s,noCollide);
        break;
    case Simulator::Algorithm::RK4Fixed:
        fast.simulateRK4Fixed(h,ts,s,noCollide);
        break;
    default:
        std::cerr<<"Unknown simulaton algorithm\n";
        break;
    }

    Simu=fast;
}


void MainWindow::drawCharts() {
    if(Simu.getResult().size()<=0) {
        std::cerr<<"No avaliable simulation result\n";
        return;
    }
//qDebug()<<__LINE__;
    const auto * result=&Simu.getResult();

//qDebug()<<__LINE__;
    QList<QVector<QPointF>> energy;
    energy.clear();
    energy.push_back(QVector<QPointF>());
    auto & energyHead=energy.first();
    energyHead.reserve(result->size());
//qDebug()<<__LINE__;

    QList<QVector<QPointF>> motions;
    for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
        motions.push_back(QVector<QPointF>());
    }
    for(auto it : motions) {
        it.clear();
        it.reserve(result->size());
    }

    for(auto it=result->cbegin();it!=result->cend();it++) {
        double yearTime=it->first/year;
        energyHead.push_back(QPointF(yearTime,
                                 Simu.calculateEnergy(it)));
        DimVector motionVal;
        Simu.calculateTotalMotion(it,motionVal);
        for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
            motions[dim].push_back(QPointF(yearTime,motionVal[dim]));
        }
    }
//qDebug()<<__LINE__;
    {
        int yPower;
        addSeriesToChart(EnergyView->chart(),
                         energy,&yPower);
        //qDebug()<<__LINE__;
        EnergyView->chart()->legend()->hide();
        EnergyView->chart()->axes(Qt::Horizontal).first()->setTitleText("Time (year)");
        EnergyView->chart()->axes(Qt::Vertical).first()->setTitleText(
                    "Energy (SI × 1E"+QString::number(yPower)+")");
        //qDebug()<<__LINE__;
        EnergyView->setRenderHint(QPainter::Antialiasing, true);
    }

    {
        int yPower;
        addSeriesToChart(MotionView->chart(),
                         motions,&yPower);
        MotionView->chart()->legend()->setAlignment(Qt::AlignmentFlag::AlignRight);
        MotionView->chart()->axes(Qt::Horizontal).first()->setTitleText("Time (year)");
        MotionView->chart()->axes(Qt::Vertical).first()->setTitleText(
                    "Motion (SI × 1E"+QString::number(yPower)+")");

        auto serieses=MotionView->chart()->series();

        for(uint16_t i=0;i<serieses.size();i++) {
            serieses[i]->setName("Dim "+QString::number(i));
        }
        MotionView->setRenderHint(QPainter::Antialiasing, true);
    }

}

void MainWindow::addSeriesToChart(QChart * chart,
                                  QList<QVector<QPointF>>& data,
                                  int * yScalePower) {
    double ymax,ymin;
    double xmin=data[0][0].x(),xmax=data[0].back().x();
    chart->removeAllSeries();
    //if(chart->series().empty()) {
        ymax=data.first()[0].y();
        ymin=ymax;
    //}
    /* else {

        const QValueAxis * xAxis=
                qobject_cast<QValueAxis*>(chart->axes(Qt::Horizontal).first());
        const QValueAxis * yAxis=
                qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).first());
        max.setX(xAxis->max());
        max.setY(xAxis->max());
        min.setX(xAxis->min());
        min.setY(yAxis->min());

        for(auto it : data) {
            min.setX(std::min(min.x(),it.x()));
            min.setY(std::min(min.y(),it.y()));

            max.setX(std::max(max.x(),it.x()));
            max.setY(std::max(max.y(),it.y()));
        }

    }*/

    for(auto it : data) {
        for(auto jt : it) {
            ymin=std::min(ymin,jt.y());
            ymax=std::max(ymax,jt.y());
        }
    }

    int yScale;
    double yRatio;
    yScale=std::floor(std::log10(
            std::max(std::abs(ymin),std::abs(ymax))));

    yRatio=std::pow(10,yScale);

    ymin/=yRatio;
    ymax/=yRatio;

    if(yScalePower!=nullptr) {
        *yScalePower=yScale;
    }

    for(auto it : data) {
        for(auto jt : it) {
            jt.setY(jt.y()/yRatio);
        }
    }

    for(auto it : data) {
        QLineSeries * series=new QLineSeries;
        series->replace(it);
        chart->addSeries(series);
    }
qDebug()<<__LINE__;
    chart->axes(Qt::Horizontal).first()->setRange(xmin,xmax);
    chart->axes(Qt::Vertical).first()->setRange(ymin,ymax);
qDebug()<<__LINE__;
}
