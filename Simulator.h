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
extern const double TimeMax; //max simulation time

class Simulator
{
public:
    Simulator();

    enum Algorithm {
        Euler,
        RK4Fixed,
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

    void clear();

    void setMass(const BodyVector &);

    const BodyVector & getMass() const;

    const std::list<Point> & getResult() const;

    double calculateKinetic(const std::_List_const_iterator<Point>) const;

    double calculatePotential(const std::_List_const_iterator<Point>) const;

    double calculateEnergy(const std::_List_const_iterator<Point>) const;

    void calculateTotalMotion(const std::_List_const_iterator<Point>,
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

private:
    std::list<Point> sol;

    BodyVector mass;


};

#endif // SIMULATOR_H
