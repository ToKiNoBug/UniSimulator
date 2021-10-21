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

    {
        auto it = pathViews.begin();
        auto dt = dimPairs.begin();
        const uint16_t colCount=std::ceil(std::sqrt(DIM_COUNT*(DIM_COUNT-1)/2));
        uint32_t chartOrder=0;
        for(uint16_t dimA=0;dimA<DIM_COUNT;dimA++) {
            for(uint16_t dimB=dimA+1;dimB<DIM_COUNT;dimB++) {
                dt->first=dimA;
                dt->second=dimB;
                *it=new QChartView(this);
                ui->dynamicLayout->addWidget(*it,
                                             chartOrder/colCount,chartOrder%colCount);
                (*it)->setChart(createEmptyChart());
                (*it)->chart()->setTitle("Path on : Dim "+QString::number(dimA)+
                                         " - Dim "+QString::number(dimB));

                (*it)->chart()->axes(Qt::Horizontal).first()->setTitleText(
                            "Dim "+QString::number(dimA));
               (*it)->chart()->axes(Qt::Vertical).first()->setTitleText(
                            "Dim "+QString::number(dimB));
                (*it)->setRenderHint(QPainter::Antialiasing, true);
                dt++;
                it++;
                chartOrder++;
            }
        }
    }

   setWindowTitle(QString::number(BODY_COUNT)+" Bodies in "
                  +QString::number(DIM_COUNT)
                  +"D Space Simulation    Made by TokiNoBug");

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
#ifdef TEST_BODY2DIM2
    BodyVector mass(Ms,2*Ms);

    Statue start;
    start.first.setValues({{-rs,rs},
                                    {0,0}});
    start.second.setValues({{-1*vs,0},
                                          {-vs,vs}});

    TimeSpan tSpan=std::make_pair(0*year,10*year);

    Time step=0.01*year;

    runSimulaton(algo,step,tSpan,start,mass);
#endif

#ifdef TEST_BODY2DIM3
    BodyVector mass(3*Ms,2*Ms);

    Statue start;
    start.first.setValues({{-rs,rs},
                                    {-0.5*rs,3*rs},
                                    {1.2*rs,5*rs}});
    start.second.setValues({{-0.9*vs,0.2*vs},
                                          {-vs,0.4*vs},
                                            {-0.3*vs,0.2*vs}});

    TimeSpan tSpan=std::make_pair(0*year,10*year);

    Time step=0.005*year;

    runSimulaton(algo,step,tSpan,start,mass);
#endif

#ifdef TEST_BODY3DIM3
    BodyVector mass(1*Ms,2*Ms,3*Ms);

    mass*=10;

    Statue start;

    Position rnd;
    rnd.setRandom();
    auto temp1=rnd*rs;
    start.first=temp1;

    rnd.setRandom();
    auto temp2=rnd*(0*vs/512);
    start.second=temp2;
    /*
    start.first.setValues({{rs,        0,     0},
                                    {0,         rs,     0},
                                    {0,         0,      rs}});
    start.second.setValues({{vs,    0,    0},
                                           {0,     vs,   0},
                                           {0,      0,      vs}});
    */
    //temp.eval();
    //start.second=temp;

    std::cout<<"starting position=\n"<<start.first<<std::endl;

    std::cout<<"start velocity=\n"<<start.second<<std::endl;

    TimeSpan tSpan=std::make_pair(0*year,5*year);

    Time step=0.001*year;

    runSimulaton(algo,step,tSpan,start,mass);
#endif
}

void MainWindow::runSimulaton(Simulator::Algorithm algo,
                                                       const double h,
                                                       const TimeSpan & ts,
                                                       const Statue & s,
                                                       const BodyVector & mv,
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


void MainWindow::drawConservativeCharts() {
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

void MainWindow::drawPathCharts() {
    if(Simu.getResult().size()<=0) {
        std::cerr<<"No avaliable simulation result\n";
        return;
    }

    std::array<std::array<QVector<QPointF>*,BODY_COUNT>,
            DIM_COUNT*(DIM_COUNT-1)/2> chartDatas;

    std::array<std::array<QLineSeries *,BODY_COUNT>,
            DIM_COUNT*(DIM_COUNT-1)/2> chartsSerieses;

    std::array<float,DIM_COUNT> posMin,posMax;//DimVector

    //double tmin=Simu.getResult().front().first,tmax=Simu.getResult().back().first;
    {
        const Point & begin=Simu.getResult().front();
        for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
            Eigen::Tensor<double,0> max,min;
            max=begin.second.first.chip(dim,0).maximum();
            min=begin.second.first.chip(dim,0).minimum();
            posMin[dim]=min(0);
            posMax[dim]=max(0);
        }
    }

    for(uint32_t chartIdx=0;chartIdx<chartDatas.size();chartIdx++) {
        for(uint32_t bodyIdx=0;bodyIdx<BODY_COUNT;bodyIdx++) {
            chartDatas[chartIdx][bodyIdx]=new QVector<QPointF>;
            chartDatas[chartIdx][bodyIdx]->resize(0);
            chartDatas[chartIdx][bodyIdx]->reserve(Simu.getResult().size());

            chartsSerieses[chartIdx][bodyIdx]=new QLineSeries;
        }
    }

    //run over the whole history
    for(auto it=Simu.getResult().cbegin();it!=Simu.getResult().cend();it++) {
        const Position & pos = it->second.first;

        //run over all dimensions to update max and min value of coordinates
        for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
            Eigen::Tensor<double,0> max,min;
            max=pos.chip(dim,0).maximum();
            min=pos.chip(dim,0).minimum();
            posMin[dim]=std::min(posMin[dim],float(min(0)));
            posMax[dim]=std::max(posMax[dim],float(max(0)));
        }

        //run over each chart
        for(uint32_t chartIdx=0;chartIdx<chartDatas.size();chartIdx++) {
            uint16_t dim_x,dim_y;
            dim_x=dimPairs[chartIdx].first;
            dim_y=dimPairs[chartIdx].second;

            //run over each body of chart
            for(uint32_t bodyIdx=0;bodyIdx<BODY_COUNT;bodyIdx++) {

                float xVal,yVal;
                xVal=pos(dim_x,bodyIdx);
                yVal=pos(dim_y,bodyIdx);

                chartDatas[chartIdx][bodyIdx]->push_back(QPointF(xVal,yVal));
            }
        }
    }

    std::array<float,DIM_COUNT> dimRatios;
    std::array<int,DIM_COUNT> dimScales;

    //run over all dimensions and process the max and min pos value
    for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
        dimScales[dim]=std::floor(std::log10(
                                      std::max(std::abs(posMin[dim]),std::abs(posMax[dim]))));
        dimRatios[dim]=std::pow(10,dimScales[dim]);

        posMin[dim]/=dimRatios[dim];
        posMax[dim]/=dimRatios[dim];
    }

    //run over all charts
    for(uint32_t chartIdx=0;chartIdx<chartDatas.size();chartIdx++) {

        uint16_t dim_x,dim_y;
        dim_x=dimPairs[chartIdx].first;
        dim_y=dimPairs[chartIdx].second;
        //run over all series (each series is a body's path)
        for(uint32_t bodyIdx=0;bodyIdx<BODY_COUNT;bodyIdx++) {
            auto curQvQpf=chartDatas[chartIdx][bodyIdx];
            //run over the whole history and scale down all values
            for(auto it=curQvQpf->begin();it!=curQvQpf->end();it++) {
                it->setX(it->x()/dimRatios[dim_x]);
                it->setY(it->y()/dimRatios[dim_y]);
            }
        }
    }

    for(uint32_t chartIdx=0;chartIdx<chartDatas.size();chartIdx++) {
        QChart * curChart=pathViews[chartIdx]->chart();
        uint16_t dim_x,dim_y;
        dim_x=dimPairs[chartIdx].first;
        dim_y=dimPairs[chartIdx].second;

        curChart->removeAllSeries();
        QValueAxis * xAxis=qobject_cast<QValueAxis *>
                (curChart->axes(Qt::Horizontal).first());
        QValueAxis * yAxis=qobject_cast<QValueAxis *>
                (curChart->axes(Qt::Vertical).first());
        xAxis->setRange(posMin[dim_x],posMax[dim_x]);
        yAxis->setRange(posMin[dim_y],posMax[dim_y]);

        for(uint32_t bodyIdx=0;bodyIdx<BODY_COUNT;bodyIdx++) {
            auto curSeries=chartsSerieses[chartIdx][bodyIdx];
            curSeries->replace(*chartDatas[chartIdx][bodyIdx]);
            curChart->addSeries(curSeries);
            curSeries->attachAxis(xAxis);
            curSeries->attachAxis(yAxis);
            delete chartDatas[chartIdx][bodyIdx];
        }
        curChart->legend()->hide();
        xAxis->setTitleText("Dim "
                            +QString::number(dim_x)
                            +" (SI × 1E"
                            +QString::number(dimScales[dim_x])+")");
        yAxis->setTitleText("Dim "
                            +QString::number(dim_y)
                            +" (SI × 1E"
                            +QString::number(dimScales[dim_y])+")");
    }



}
