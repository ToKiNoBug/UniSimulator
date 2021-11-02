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

#include "MainWindow.h"
#include "ui_MainWindow.h"

const QColor DefaultColors[]={
                       Qt::GlobalColor::blue,
                       Qt::GlobalColor::green,
                       Qt::GlobalColor::yellow,
                       Qt::GlobalColor::magenta,
                       Qt::GlobalColor::red,
                       Qt::GlobalColor::gray,
                       Qt::GlobalColor::cyan,
                       Qt::GlobalColor::darkGreen};

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow) {
    ui->setupUi(this);

    buildMotionEnergyUI();

    buildPathUI();

    buildParameterUI();

   setWindowTitle(QString::number(BODY_COUNT)+" Bodies in "
                  +QString::number(DIM_COUNT)
                  +"D Space Simulation    Made by TokiNoBug");
    isRunning=false;

    connect(ui->actionSave_parameters,&QAction::triggered,
            this,&MainWindow::saveParameters);
    connect(ui->actionLoad_parameters,&QAction::triggered,
            this,&MainWindow::loadParameters);
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

    Time step=0.001*year;

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

#endif

#ifdef TEST_BODY3DIM3
    BodyVector mass(1*Ms,2*Ms,3*Ms);

    mass*=1;

    Statue start;

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

    std::cout<<"starting position=\n"<<start.first/rs<<std::endl;

    std::cout<<"starting velocity=\n"<<start.second/vs<<std::endl;

    TimeSpan tSpan=std::make_pair(0*year,100*year);

    Time step=1e-4*year;
    setParamaters(mass,start,tSpan,step,algo);
    on_BtnRunSimulation_clicked();

#endif
}

void MainWindow::runSimulaton(Simulator::Algorithm algo,
                                                       const double h,
                                                       const TimeSpan & ts,
                                                       const Statue & s,
                                                       const BodyVector & mv,
                                                       bool * noCollide) {
    isRunning=true;
    Simulator fast=Simu;

    fast.setMass(mv);

    switch (algo) {
    case Simulator::Algorithm::Euler:
        fast.simulateEuler(h,ts,s,noCollide);
        break;
    case Simulator::Algorithm::RK4Fixed:
        fast.simulateRK4Fixed(h,ts,s,noCollide);
        break;
    case Simulator::Algorithm::RK4Var1:
        fast.simulateRK4Var1(h,ts,s,noCollide);
        break;
    default:
        std::cerr<<"Unknown simulaton algorithm\n";
        return;
    }
    std::cerr<<"Simulation finished with "<<fast.getResult().size()<<" dots\n";
    //

    uint64_t dotCount;
#ifdef SIMU_INTERPOLT
    double showBeg,showEnd;
    showBeg=ts.first;
    showEnd=std::min(ts.second,fast.getResult().back().first);
    dotCount=std::max(512,int((showEnd-showBeg)*32/year));

    Simulator::deval(&fast,&Simu,
                     Eigen::ArrayXd::LinSpaced(dotCount,showBeg,showEnd));
#else
    Simu=fast;
    dotCount=fast.getResult().size();
#endif
    std::cerr<<"curve displayed with "<<dotCount<<" dots\n";
    ui->timeSlider->setMinimum(0);
    ui->timeSlider->setMaximum(dotCount-1);

    isRunning=false;
}

void MainWindow::on_timeSlider_valueChanged(int value) {

    if(isRunning) return;

    //std::cerr<<"value="<<value<<std::endl;
    for(auto it : pathViews) {
        moveScatterIndex(it->chart(),value);
    }

    //return;
    //qDebug()<<__LINE__;
    moveScatterIndex(MotionView->chart(),value);
    QPointF p=moveScatterIndex(EnergyView->chart(),value);
    ui->displayTime->setText("Time : "+QString::number(p.x())+" year(s)");
}


void MainWindow::on_timeSlider_sliderPressed() {

    QChart::AnimationOptions op;

    op.setFlag(QChart::AnimationOption::GridAxisAnimations,false);
    op.setFlag(QChart::AnimationOption::SeriesAnimations,false);

    for(auto it : pathViews) {
        it->chart()->setAnimationOptions(op);
    }
    EnergyView->chart()->setAnimationOptions(op);
    MotionView->chart()->setAnimationOptions(op);
}


void MainWindow::on_timeSlider_sliderReleased() {
    QChart::AnimationOptions op;

    op.setFlag(QChart::AnimationOption::GridAxisAnimations,false);
    op.setFlag(QChart::AnimationOption::SeriesAnimations,true);

    for(auto it : pathViews) {
        it->chart()->setAnimationOptions(op);
    }
    EnergyView->chart()->setAnimationOptions(op);
    MotionView->chart()->setAnimationOptions(op);

}

void MainWindow::on_cbShowCurves_stateChanged(int arg1) {
    std::list<QAbstractSeries*> series;

    series.resize(0);

    for(auto it : pathViews) {
        //=it->chart()->series();
        QList<QAbstractSeries*> temp=it->chart()->series();
        for(auto jt : temp) {
            if(jt->type()==QAbstractSeries::SeriesType::SeriesTypeLine) {
                series.push_back(jt);
            }
        }
    }

    for(auto jt : EnergyView->chart()->series()) {
        if(jt->type()==QAbstractSeries::SeriesType::SeriesTypeLine) {
            series.push_back(jt);
        }
    }

    for(auto jt : MotionView->chart()->series()) {
        if(jt->type()==QAbstractSeries::SeriesType::SeriesTypeLine) {
            series.push_back(jt);
        }
    }

    for(auto jt : series) {
            if(arg1) {
                jt->show();
            } else {
                jt->hide();
            }
    }

}

void MainWindow::on_cbShowDots_stateChanged(int arg1) {
    std::list<QAbstractSeries*> series;

    series.resize(0);

    for(auto it : pathViews) {
        //=it->chart()->series();
        QList<QAbstractSeries*> temp=it->chart()->series();
        for(auto jt : temp) {
            if(jt->type()==QAbstractSeries::SeriesType::SeriesTypeScatter) {
                series.push_back(jt);
            }
        }
    }

    for(auto jt : EnergyView->chart()->series()) {
        if(jt->type()==QAbstractSeries::SeriesType::SeriesTypeScatter) {
            series.push_back(jt);
        }
    }

    for(auto jt : MotionView->chart()->series()) {
        if(jt->type()==QAbstractSeries::SeriesType::SeriesTypeScatter) {
            series.push_back(jt);
        }
    }

    for(auto jt : series) {
            if(arg1) {
                jt->show();
            } else {
                jt->hide();
            }
    }
}

#include "WidgetCodes.cpp"
