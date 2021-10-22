#include "MainWindow.h"

#include <QApplication>
#include <ctime>

#include "tests.h"

int main(int argc, char *argv[])
{

    dispConstants();

    QApplication a(argc, argv);
    MainWindow w;
    w.show();
#ifdef TEST_BODY2DIM2
    w.runSimulaton(Simulator::Algorithm::RK4Var1);
    w.drawConservativeCharts();
    w.drawPathCharts();
#endif

#ifdef TEST_BODY2DIM3
    w.runSimulaton(Simulator::Algorithm::RK4Fixed);
    w.drawConservativeCharts();
    w.drawPathCharts();
#endif

#ifdef TEST_BODY3DIM3
    w.runSimulaton(Simulator::Algorithm::RK4Var1);
    w.drawConservativeCharts();
    w.drawPathCharts();
#endif

    return a.exec();
}
