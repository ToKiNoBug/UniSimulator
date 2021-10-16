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

    w.runSimulaton(Simulator::Algorithm::RK4Fixed);
    w.drawCharts();

    return a.exec();
}
