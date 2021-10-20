#ifndef DEFINES_H
#define DEFINES_H


//#define EIGEN_NO_DEBUG

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#define DIM_COUNT 3
#define BODY_COUNT 2

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

typedef Eigen::Array<double,DIM_COUNT,1> DimVector ;

typedef Eigen::Array<double,BODY_COUNT,BODY_COUNT> DistanceMat ;

typedef std::pair<Position,Velocity> Statue ;

typedef  std::pair<Velocity,Acceleration> Derivative ;

typedef std::pair<Time,Time> TimeSpan ;

typedef  std::pair<Time,Statue> Point ;


#if DIM_COUNT <=0
Error : DIM_COUNT_should_be_a_positive_integer
#endif

#if BODY_COUNT <2
Error : BODY_COUNT_should_be_a_positive_integer_not_less_than_2
#endif


#endif // DEFINES_H
