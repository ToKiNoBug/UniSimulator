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

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <tuple>
#include <fstream>
#include <string>

#include "defines.h"

extern const double G; //gravity constant
extern const double Ms; //standard mass (solar mass)
extern const double omega_s; //standard angle speed
extern const double rho; //solar density
extern const double year; //seconds of a year
extern const double rs; //standard distance (1AU)
extern const double vs; //standard speed
extern const double as; //standard accelerate

class Simulator
{
public:
    Simulator();

    enum Algorithm {
        Euler,
        RK4Fixed,
        RK4Var1,
        //ODE45
    };
    static const std::string paraSuffix;
    static const std::string dataSuffix;
public:
    void simulateEuler(const double,
                       TimeSpan,
                       Statue,
                       bool * noCollide=nullptr);

    void simulateRK4Fixed(const double,
                          TimeSpan,
                          Statue,
                          bool * noCollide=nullptr);

    void simulateRK4Var1(double,
                         TimeSpan,
                         Statue,
                         bool * noCollide=nullptr);

    void clear();

    void setMass(const BodyVector &);

    const BodyVector & getMass() const;

    const std::list<Point> & getResult() const;

    double calculateKinetic(const Statue & it) const;

    double calculatePotential(const Statue & it) const;

    double calculateEnergy(const Statue & it) const;

    void calculateTotalMotion(const Statue & it,
                                                DimVector &) const;

    void savePath(const char * fileName);

public:
    static void calculateSafeDistance(const BodyVector &,
                                                            DistanceMat & dest);

    static void calculateGM(const BodyVector &,
                                                            Interaction &);

    static bool calculateDiff(const Position & y,
                             const Interaction &,
                             const DistanceMat &,
                             Acceleration & dy);
    //To avoid useless deep copying, dy.first will not be used.

    static bool RK4(const Time h,
                    const Statue & y,
                    const Interaction &,
                    const DistanceMat& ,
                    Statue & y_n1);

    static bool isErrorTolerantable(const Statue & y_h,
                                    const Statue & y_h_2,
                                    double errorRatio=1e-8);

    static void deval(const Simulator * source,
            Simulator * dest,
                          const Eigen::ArrayXd & timeQueried);
    static void motionAlign(const BodyVector & mass, Velocity & velocity);
    static void positonAlign(Position &);
    static void saveParameters(const char * fileName,
                               const BodyVector & mass,
                               const Statue & y0,
                               const TimeSpan ts,
                               const double step);
    static bool loadParameters(const char * fileName,
                               BodyVector & mass,
                               Statue & y0,
                               TimeSpan & ts,
                               double & step);
    void saveAsData(const char * fileName) const;
    bool loadFromData(const char * fileName);

#ifdef BODY3_DIM3
    static void positionShrink(Position &);
#endif
private:
    std::list<Point> sol;

    BodyVector mass;


};

#endif // SIMULATOR_H
