#include "MainWindow.h"

#include <QApplication>
#include <ctime>

#if (BODY_COUNT==2) && (DIM_COUNT==2)
#define TEST
#define TEST_BODY2DIM2
#endif

#ifdef TEST
void testDerivative();
void testPerformance();
void testSimulation();
#endif

int main(int argc, char *argv[])
{

    //testSimulation();

    QApplication a(argc, argv);
    MainWindow w;
    w.show();



    return a.exec();
}


#ifdef TEST_BODY2DIM2
void testPerformance() {
std::cerr<<"test started";

MassVector M(Ms,1*Ms);

Statue y,dy;
y.first.setValues({{-rs,rs},
                            {0,0}});
y.second.setValues({{0,0},
                                {-vs,vs}});

DistanceMat safeMat;
Simulator::calculateSafeDistance(M,safeMat);

std::cerr<<"safe distance=\n";
std::cerr<<safeMat<<std::endl;

Interaction GM;
Simulator::calculateGM(M,GM);

clock_t curTime=std::clock();

int repeatTime=100000;

double h=1e-2*year;

for(int i=0;i<repeatTime;i++) {
    if(i%1024==0) {

        std::cerr<<i<<std::endl;

        std::cerr<<"position=\n";
        std::cerr<<y.first<<std::endl;

        std::cerr<<"velocity=\n";
        std::cerr<<y.second<<std::endl;
    }
//bool isOk=Simulator::calculateDiff(y.first,GM,safeMat,dy.second);
    bool isOk=Simulator::RK4(h,y,GM,safeMat,y);
if(!isOk) {
    std::cerr<<"stars will collide\n";
    break;
}
/*
y.first+=h*y.second;
y.second+=h*dy.second;
*/
}

std::cerr<<double(clock()-curTime)*1e6/repeatTime/CLOCKS_PER_SEC
             <<"*10^-6 s"<<std::endl;
/*
std::cerr<<"acceleration=\n";
std::cerr<<dy.second<<std::endl;
*/

//std::cerr<<(isOk?"stars won't collide\n":"stars will collide\n");

exit(0);
return;
}

void testDerivative() {
    std::cout<<"testDerivative\n";

    MassVector M(Ms,Ms);

    Statue y;
    y.first.setValues({{-rs,rs},
                                {0,0}});
    y.second.setValues({{0,0},
                                    {-vs,vs}});

    Acceleration acc;

    DistanceMat safeMat;
    Simulator::calculateSafeDistance(M,safeMat);


    Interaction GM;
    Simulator::calculateGM(M,GM);

    std::cerr<<"\nposition:\n";
    std::cerr<<y.first;
    std::cerr<<"\nvelocity:\n";
    std::cerr<<y.second;
    Simulator::calculateDiff(y.first,GM,safeMat,acc);

    std::cerr<<"\nacceleration:\n";
    std::cerr<<acc;

    std::cout<<"\nfinished\n";
    exit(0);
    return;
}

void testSimulation() {
    std::cout<<"testSimulation\n";

    MassVector M(Ms,Ms);

    Statue y;
    y.first.setValues({{-rs,rs},
                                {0,0}});
    y.second.setValues({{0,0},
                                    {-vs,vs}});

    Acceleration acc;

    DistanceMat safeMat;
    Simulator::calculateSafeDistance(M,safeMat);


    Interaction GM;
    Simulator::calculateGM(M,GM);

    Time maxTime=100*year;
    Time stepL=1e-2*year;

    uint64_t loopCount=maxTime/stepL;

    for(uint64_t i=0;i<loopCount;i++) {
        std::cerr<<"\nposition:\n";
        std::cerr<<y.first;
        std::cerr<<"\nvelocity:\n";
        std::cerr<<y.second;
        std::cerr<<"\nacceleration:\n";
        Simulator::calculateDiff(y.first,GM,safeMat,acc);
        std::cerr<<acc;


        bool isOk=Simulator::RK4(stepL,y,GM,safeMat,y);
        if(!isOk) {
            std::cerr<<"\nStars will collide\n";
            break;
        }
    }
    std::cout<<"\nfinished\n";
    exit(0);
    return;
}
#endif
