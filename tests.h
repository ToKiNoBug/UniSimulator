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

#ifndef TESTS_H
#define TESTS_H

#include "Simulator.h"

#if (BODY_COUNT==2) && (DIM_COUNT==2)
#define TEST_BODY2DIM2
#endif

#if (BODY_COUNT==2) && (DIM_COUNT==3)
#define TEST_BODY2DIM3
#endif

#if (BODY_COUNT==3) && (DIM_COUNT==3)
#define TEST_BODY3DIM3
#endif

void dispConstants();

#ifdef TEST_BODY2DIM2
void testDerivative();
void testPerformance();
void testSimulation();

void testEuler();
void testRK4Fixed();

#endif

#endif // TESTS_H
