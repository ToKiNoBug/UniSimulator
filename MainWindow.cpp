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
    MassVector mass(Ms,2*Ms);

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
    QList<QVector<QPointF>> energy(1);
    auto & energyHead=energy.first();
    energyHead.resize(0);
    energyHead.reserve(result->size());
//qDebug()<<__LINE__;

    QList<QVector<QPointF>> motions(DIM_COUNT);

    for(auto it : motions) {
        it.clear();
        it.resize(0);
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

        qDebug()<<__LINE__;
        addSeriesToChart(MotionView->chart(),
                         motions,&yPower);
        qDebug()<<__LINE__;
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
    const double xmin=data[0].first().x(),xmax=data[0].back().x();
    double ymax,ymin;
    chart->removeAllSeries();
    //if(chart->series().empty()) {
        ymax=data[0][0].y();
        ymin=ymax;
    //}
    for(auto it=data.cbegin();it!=data.cend();it++) {
        for(auto jt=it->cbegin();jt!=it->cend();jt++) {
            ymin=std::min(ymin,jt->y());
            ymax=std::max(ymax,jt->y());
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

    for(auto it=data.begin();it!=data.end();it++) {
        for(auto jt=it->begin();jt!=it->end();jt++) {
            jt->setY(jt->y()/yRatio);
        }
    }

    QValueAxis * xAxis=qobject_cast<QValueAxis *>(chart->axes(Qt::Horizontal).first());
    QValueAxis * yAxis=qobject_cast<QValueAxis *>(chart->axes(Qt::Vertical).first());

    qDebug()<<__LINE__;
        xAxis->setRange(xmin,xmax);
        yAxis->setRange(ymin,ymax);

    qDebug()<<__LINE__;
    for(uint16_t d=0;d<data.size();d++) {
        QLineSeries * series=new QLineSeries;
        series->replace(data[d]);
        chart->addSeries(series);

        series->attachAxis(xAxis);
        series->attachAxis(yAxis);
    }

}
