#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <cmath>
#include <list>
#include <vector>

//#define EIGEN_NO_DEBUG

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#define DIM_COUNT 3
#define BODY_COUNT 2

extern const double G;
extern const double Ms;
extern const double omega_s;
extern const double rho;
extern const double year;
extern const double rs;
extern const double vs;
extern const double as;
extern const double TimeMax;

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

typedef std::pair<Time,Statue> Point ;

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

    static bool calculateDiff(const Statue & y,
                             const Interaction &,
                             const DistanceMat &,
                             Acceleration & dy);
    //To avoid useless deep copying, dy.first will not be used.

    static bool RK4(const double h,
                    const Statue & y,
                    const Interaction &,
                    const DistanceMat& ,
                    Statue & y_n1);

private:
    std::list<Point> sol;


};

#endif // SIMULATOR_H
