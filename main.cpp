/*
 Copyright © 2021  TokiNoBug
This file is part of UniSimulator.

    UniSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    UniSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with UniSimulator.  If not, see <https://www.gnu.org/licenses/>.

    Contact with me:
    github:https://github.com/ToKiNoBug
    bilibili:https://space.bilibili.com/351429231
*/

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
#endif

#ifdef TEST_BODY2DIM3
    w.runSimulaton(Simulator::Algorithm::RK4Fixed);
#endif

#ifdef TEST_BODY3DIM3
    w.runSimulaton(Simulator::Algorithm::RK4Var1);
#endif

    return a.exec();
}
