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
                                 Simu.calculateEnergy(it->second)));
        DimVector motionVal;
        Simu.calculateTotalMotion(it->second,motionVal);
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

        //qDebug()<<__LINE__;
        addSeriesToChart(MotionView->chart(),
                         motions,&yPower);
        //qDebug()<<__LINE__;
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

    //qDebug()<<__LINE__;
        xAxis->setRange(xmin,xmax);
        yAxis->setRange(ymin,ymax);

    //qDebug()<<__LINE__;
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
        //std::cerr<<"Line="<<__LINE__<<"chartIdx="<<chartIdx<<std::endl;
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

void MainWindow::buildParameterUI() {

    connect(ui->inputBeginTime,&QLineEdit::textChanged,
            this,&MainWindow::onUserInput);
    connect(ui->inputEndTime,&QLineEdit::textChanged,
            this,&MainWindow::onUserInput);
    connect(ui->inputStep,&QLineEdit::textChanged,
            this,&MainWindow::onUserInput);

    ui->massBox->setTitle("Mass ( Ms="+QString::number(Ms)+"kg )");
    for(uint32_t body=0;body<BODY_COUNT;body++) {
        QLabel * curLabel=new QLabel;
        curLabel->setText("M"+QString::number(body)+"=");
        ui->massArea->addWidget(curLabel,0,2*body);
        QLineEdit * curWidget=new QLineEdit;
        massWidgets[body]=curWidget;
        ui->massArea->addWidget(curWidget,0,2*body+1);
        connect(curWidget,&QLineEdit::textChanged,
                this,&MainWindow::onUserInput);
    }

    ui->positionBox->setTitle("Position ( rs="+QString::number(rs)+"m )");
    for(uint32_t dim=0;dim<DIM_COUNT;dim++) {
        for(uint32_t body=0;body<BODY_COUNT;body++) {
            QLineEdit * curWidget=new QLineEdit;
            ui->positionArea->addWidget(curWidget,dim,body);
            positionWidgets[dim][body]=curWidget;
            connect(curWidget,&QLineEdit::textChanged,
                    this,&MainWindow::onUserInput);
        }
    }

    ui->velocityBox->setTitle("Velocity ( vs="+QString::number(vs)+"m/s )");
    for(uint32_t dim=0;dim<DIM_COUNT;dim++) {
        for(uint32_t body=0;body<BODY_COUNT;body++) {
            QLineEdit * curWidget=new QLineEdit;
            ui->velocityArea->addWidget(curWidget,dim,body);
            velocityWidgets[dim][body]=curWidget;
            connect(curWidget,&QLineEdit::textChanged,
                    this,&MainWindow::onUserInput);
        }
    }

    ui->timeBox->setTitle("Time ( year="+QString::number(year)+"s )");

    ui->selectAlgorithm->addItem("Euler",Simulator::Algorithm::Euler);
    ui->selectAlgorithm->addItem("RK4Fixed",Simulator::Algorithm::RK4Fixed);
    ui->selectAlgorithm->addItem("RK4Var1",Simulator::Algorithm::RK4Var1);
}

void MainWindow::setParamaters(const BodyVector & mass,
                                const Statue & y0,
                               TimeSpan ts,
                               double step,
                               Simulator::Algorithm algo) {
    for(uint32_t body=0;body<BODY_COUNT;body++) {
        massWidgets[body]->setText(QString::number(mass[body]/Ms));
    }

    for(uint32_t dim=0;dim<DIM_COUNT;dim++) {
        for(uint32_t body=0;body<BODY_COUNT;body++) {
            positionWidgets[dim][body]->setText(QString::number(y0.first(dim,body)/rs));
        }
    }

    for(uint32_t dim=0;dim<DIM_COUNT;dim++) {
        for(uint32_t body=0;body<BODY_COUNT;body++) {
            velocityWidgets[dim][body]->setText(QString::number(y0.second(dim,body)/vs));
        }
    }

    ui->inputBeginTime->setText(QString::number(ts.first/year));
    ui->inputEndTime->setText(QString::number(ts.second/year));
    ui->inputStep->setText(QString::number(step/year));

    uint16_t selectedIdx=0;
    switch (algo) {
    case Simulator::Algorithm::Euler:
        selectedIdx=0;
        break;
    case Simulator::Algorithm::RK4Fixed:
        selectedIdx=1;
        break;
    case Simulator::Algorithm::RK4Var1:
        selectedIdx=2;
        break;
    default:
        selectedIdx=2;
        break;
    }
    ui->selectAlgorithm->setCurrentIndex(selectedIdx);
}


bool MainWindow::grabParameters(BodyVector &mass,
                                Statue &y0,
                                TimeSpan &ts,
                                double &step,
                                Simulator::Algorithm & algo) {
    bool ok=true;
    for(uint32_t body=0;body<BODY_COUNT;body++) {
        mass[body]=Ms*this->massWidgets[body]->text().toDouble(&ok);
        if(!ok) return false;
        for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
            y0.first(dim,body)=rs*this->positionWidgets[dim][body]->text().toDouble(&ok);
            if(!ok) return false;

            y0.second(dim,body)=vs*this->velocityWidgets[dim][body]->text().toDouble(&ok);
            if(!ok) return false;
        }
    }

    ts.first=year*ui->inputBeginTime->text().toDouble(&ok);
    if(!ok) return false;
    ts.second=year*ui->inputEndTime->text().toDouble(&ok);
    if(!ok) return false;
    step=year*ui->inputStep->text().toDouble(&ok);
    if(!ok) return false;
    algo=Simulator::Algorithm(ui->selectAlgorithm->currentData().toInt(&ok));

    return ok;
}

void MainWindow::onUserInput() {
    double step;
    TimeSpan ts;
    Statue y0;
    BodyVector mass;
    Simulator::Algorithm algo;

    bool ok=grabParameters(mass,y0,ts,step,algo);
    ui->BtnRunSimulation->setEnabled(ok);
}


void MainWindow::on_BtnRunSimulation_clicked() {
    double step;
    TimeSpan ts;
    Statue y0;
    BodyVector mass;
    Simulator::Algorithm algo;

    bool ok=grabParameters(mass,y0,ts,step,algo);
    if(!ok) return;
    std::cerr<<"mass=\n"<<mass/Ms<<std::endl;
    std::cerr<<"position=\n"<<y0.first/rs<<std::endl;
    std::cerr<<"velocity=\n"<<y0.second/vs<<std::endl;
    std::cerr<<"time span=["<<ts.first/year<<" , "<<ts.second/year<<"]\n";
    std::cerr<<"time step="<<step/year<<std::endl;
    std::cerr<<"algo : ";
    switch (algo) {
    case Simulator::Algorithm::Euler:
        std::cerr<<"Euler\n";
        break;
    case Simulator::Algorithm::RK4Fixed:
        std::cerr<<"RK4Fixed\n";
        break;
    case Simulator::Algorithm::RK4Var1:
        std::cerr<<"RK4Var1\n";
        break;
    default:
        std::cerr<<"Unknown algorithm, use RK4Var1 as default\n";
        algo=Simulator::Algorithm::RK4Var1;
        break;
    }
    bool noCollide=true;
    runSimulaton(algo,step,ts,y0,mass,&noCollide);
    drawConservativeCharts();
    drawPathCharts();

    QString title=QString::number(BODY_COUNT)+" Bodies in "
            +QString::number(DIM_COUNT)
            +"D Space Simulation    Made by TokiNoBug";
    if(!noCollide) {
        title+=" | Stars will Collide at "
                +QString::number(Simu.getResult().back().first/year)
                +" year(s)";
    }
    this->setWindowTitle(title);

    ui->tabWidget->setCurrentIndex(1);
}


void MainWindow::on_BtnClear_clicked() {
for(auto it : massWidgets) {
    it->clear();
}
for(auto & it : positionWidgets) {
    for(auto jt : it) {
        jt->clear();
    }
}
for(auto & it : velocityWidgets) {
    for(auto jt : it) {
        jt->clear();
    }
}

ui->inputBeginTime->setText(QString::number(0));
ui->inputEndTime->setText(QString::number(100));
ui->inputStep->setText(QString::number(0.0001));

ui->selectAlgorithm->setCurrentIndex(2);

}


void MainWindow::on_BtnRandom_clicked() {
    BodyVector mass;
    Statue start;

    mass.setRandom();

    auto abs=5*mass.abs();
    abs.eval();
    mass=abs*Ms;


    Position rnd;
    rnd.setRandom();
    for(uint32_t i=0;i<rnd.size();i++){
        rnd(i)=2*(rnd(i)-0.5);
    }
    auto temp1=rnd*(10*rs);
    start.first=temp1;

    rnd.setRandom();
    for(uint32_t i=0;i<rnd.size();i++){
        rnd(i)=2*(rnd(i)-0.5);
    }



    auto temp2=(rnd)*(vs);
    start.second=temp2;
    TimeSpan tSpan=std::make_pair(0*year,100*year);

    Time step=1e-4*year;

    Simulator::Algorithm algo=Simulator::Algorithm::RK4Var1;


    setParamaters(mass,start,tSpan,step,algo);
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
    QWidget::keyPressEvent(event);
    if(event->key()==Qt::Key::Key_F5) {
        QPixmap pixmap=this->grab();

        QString path;
        path=QFileDialog::getSaveFileName(this,
                                          "Save screenshot to image",
                                          "",
                                          "*.png"
                                          );
        if(path.isEmpty()) {
            return;
        }
        pixmap.save(path);
    }
}

void MainWindow::saveParameters() {
    Statue start;
    TimeSpan ts;
    BodyVector mass;
    double step;
    Simulator::Algorithm algo;
    bool ok=
            grabParameters(mass,start,ts,step,algo);

    if(!ok) {
        return;
    }

    QString fileName=QFileDialog::getSaveFileName(this,
                                                  "Save parameters to file",
                                                  "",
                                                  "*"+QString::fromStdString(Simulator::paraSuffix));
    if(fileName.isEmpty()) {
        return;
    }

    Simulator::saveParameters(fileName.toLocal8Bit(),mass,start,ts,step);

}

void MainWindow::loadParameters() {

    Statue start;
    TimeSpan ts;
    BodyVector mass;
    double step;

    QString fileName=QFileDialog::getOpenFileName(this,
                                                  "Load parameters from file",
                                                  "",
                                                  "*"+QString::fromStdString(Simulator::paraSuffix));

    if(fileName.isEmpty()) {
        return;
    }

    bool ok=Simulator::loadParameters(fileName.toLocal8Bit(),mass,start,ts,step);

    if(!ok) {
        std::cerr<<"Failed to load parameters\n";
        return;
    }

    setParamaters(mass,start,ts,step,Simulator::Algorithm::RK4Var1);
}
#endif
