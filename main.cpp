#include "MainWindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    std::cerr<<"test started";

    MassVector M(2e30,2e30);

    return a.exec();
}
