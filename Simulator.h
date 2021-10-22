#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <tuple>


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

    double calculateKinetic(const Point*) const;

    double calculatePotential(const Point* it) const;

    double calculateEnergy(const Point* it) const;

    void calculateTotalMotion(const Point* it,
                                                DimVector &) const;


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
private:
    std::list<Point> sol;

    BodyVector mass;


};

#endif // SIMULATOR_H
