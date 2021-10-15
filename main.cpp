#include "MainWindow.h"

#include <QApplication>
#include <ctime>

void test();

int main(int argc, char *argv[])
{
    test();

    QApplication a(argc, argv);
    MainWindow w;
    w.show();



    return a.exec();
}

void test() {
std::cerr<<"test started";

MassVector M(Ms,1*Ms);

Statue y,dy;
y.first.setValues({{-rs,rs},{-rs,rs},{-rs,rs}});
y.second.setValues({{-vs,vs},{-vs,vs},{-vs,vs}});

std::cerr<<"position=\n";
std::cerr<<y.first<<std::endl;

std::cerr<<"velocity=\n";
std::cerr<<y.second<<std::endl;

DistanceMat safeMat;
Simulator::calculateSafeDistance(M,safeMat);

std::cerr<<"safe distance=\n";
std::cerr<<safeMat<<std::endl;

Interaction GM;
Simulator::calculateGM(M,GM);

clock_t curTime=std::clock();

int repeatTime=100000;

double h=1e-3*year;

for(int i=0;i<repeatTime;i++) {
    if(i%1024==0)
std::cerr<<i<<std::endl;
bool isOk=Simulator::calculateDiff(y,GM,safeMat,dy.second);
if(!isOk) {
    std::cerr<<"stars will collide\n";
    break;
}
y.first+=h*y.second;
y.second+=h*dy.second;
}

std::cerr<<double(clock()-curTime)*1e6/repeatTime/CLOCKS_PER_SEC
             <<"*10^-6 s"<<std::endl;
std::cerr<<"acceleration=\n";
std::cerr<<dy.second<<std::endl;

//std::cerr<<(isOk?"stars won't collide\n":"stars will collide\n");

exit(0);
return;
}
