#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <tuple>

//#define EIGEN_NO_DEBUG

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#define DIM_COUNT 2
#define BODY_COUNT 2

extern const double G; //gravity constant
extern const double Ms; //standard mass (solar mass)
extern const double omega_s; //standard angle speed
extern const double rho; //solar density
extern const double year; //seconds of a year
extern const double rs; //standard distance (1AU)
extern const double vs; //standard speed
extern const double as; //standard accelerate
extern const double TimeMax; //max simulation time

#ifdef IS_POINT_TUPLE
extern const uint8_t tupleTimeIndex,
                                tuplePositionIndex,
                                tupleVelocityIndex,
                                tupleAccelerlationIndex;
#endif

//typedef Eigen::TensorFixedSize<double,Eigen::Sizes<DimCount,BodyCount>> Position;

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>> Position;

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>> Velocity;

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>> Acceleration;

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT,BODY_COUNT>> Interaction;

typedef double Time ;

typedef Eigen::Array<double,BODY_COUNT,1> MassVector ;

typedef Eigen::Array<double,BODY_COUNT,BODY_COUNT> DistanceMat ;

typedef std::pair<Position,Velocity> Statue ;

typedef  std::pair<Velocity,Acceleration> Derivative ;

#ifdef IS_POINT_TUPLE
typedef std::tuple<Time,Position,Velocity,Acceleration> Point;
#else
typedef  std::pair<Time,Statue> Point ;
#endif

class Simulator
{
public:
    Simulator();

    void simulate(const Statue &);

    const std::list<Point> & getResult() const;

    static void calculateSafeDistance(const MassVector &,
                                                            DistanceMat & dest);

    static void calculateGM(const MassVector &,
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


};

#endif // SIMULATOR_H
