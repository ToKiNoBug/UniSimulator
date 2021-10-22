#ifndef WIDGETCODES_CPP
#define WIDGETCODES_CPP

#include "MainWindow.h"
#include "ui_MainWindow.h"

void MainWindow::buildMotionEnergyUI() {
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
}

void MainWindow::buildPathUI() {
    auto it = pathViews.begin();
    auto dt = dimPairs.begin();
    const uint16_t colCount=std::ceil(std::sqrt(pathChartCount));
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
                                 Simu.calculateEnergy(&*it)));
        DimVector motionVal;
        Simu.calculateTotalMotion(&*it,motionVal);
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

    createScatters(MotionView->chart());
    createScatters(EnergyView->chart());

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
            pathChartCount> chartDatas;

    std::array<std::array<QLineSeries *,BODY_COUNT>,
            pathChartCount> chartsSerieses;

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
        createScatters(curChart);
        std::cerr<<"Line="<<__LINE__<<"chartIdx="<<chartIdx<<std::endl;
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

void MainWindow::createScatters(QChart * chart) {
    uint32_t dataCount=chart->series().size();
    if(dataCount<=0)
        return;

    std::vector<QLineSeries*> series;
    series.resize(0);

    for(auto it : chart->series()) {
        series.push_back(qobject_cast<QLineSeries*>(it));
    }

    for(uint32_t dataIdx=0;dataIdx<dataCount;dataIdx++) {
        QScatterSeries * ss=new QScatterSeries;
        ss->setColor(series[dataIdx]->pen().color());
        ss->clear();
        ss->append(series[dataIdx]->at(0));
        chart->addSeries(ss);
        ss->attachAxis(chart->axes(Qt::Horizontal).front());
        ss->attachAxis(chart->axes(Qt::Vertical).front());
    }
}

QPointF MainWindow::moveScatterIndex(QChart * chart, int index) {
    int dataCount=chart->series().size();
    if((dataCount<=0)||(dataCount%2!=0))
        return QPointF();
    dataCount/=2;
    QLineSeries* ls;
    QScatterSeries* ss;
    QPointF result;
    for(int dataIdx=0;dataIdx<dataCount;dataIdx++) {
        ls=qobject_cast<QLineSeries*>(chart->series()[dataIdx+0*dataCount]);
        ss=qobject_cast<QScatterSeries*>(chart->series()[dataIdx+1*dataCount]);
        ss->clear();
        ss->append(result=ls->at(index));
        //ss->append(ls->at(std::min(index,ls->count()-1)));
    }
    return result;
}

#endif
