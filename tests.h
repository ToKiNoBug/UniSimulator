#ifndef TESTS_H
#define TESTS_H

#include "Simulator.h"

#if (BODY_COUNT==2) && (DIM_COUNT==2)
#define TEST_BODY2DIM2
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
