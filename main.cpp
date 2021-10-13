#include "MainWindow.h"

#include <QApplication>

void test();

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    test();

    return a.exec();
}

void test() {
std::cerr<<"test started";

MassVector M(Ms,1*Ms);

Statue y,dy;
y.first.setValues({{-rs,rs}});
y.second.setValues({{-vs,vs}});

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

bool isOk=Simulator::calculateDiff(y,GM,safeMat,dy);

std::cerr<<"acceleration=\n";
std::cerr<<dy.second<<std::endl;

std::cerr<<(isOk?"stars won't collide\n":"stars will collide\n");
}
