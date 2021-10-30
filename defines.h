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

#ifndef DEFINES_H
#define DEFINES_H


//#define EIGEN_NO_DEBUG

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#define DIM_COUNT 4
#define BODY_COUNT 5

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>> Position;

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>> Velocity;

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT>> Acceleration;

typedef Eigen::TensorFixedSize<double,
                        Eigen::Sizes<DIM_COUNT,BODY_COUNT,BODY_COUNT>> Interaction;

typedef double Time ;

typedef Eigen::Array<double,BODY_COUNT,1> BodyVector ;

typedef Eigen::Array<double,DIM_COUNT,1> DimVector ;

typedef Eigen::Array<double,BODY_COUNT,BODY_COUNT> DistanceMat ;

typedef std::pair<Position,Velocity> Statue ;

typedef  std::pair<Velocity,Acceleration> Derivative ;

typedef std::pair<Time,Time> TimeSpan ;

typedef  std::pair<Time,Statue> Point ;


#if (BODY_COUNT==2) && (DIM_COUNT==2)
#define TEST_BODY2DIM2
#endif

#if (BODY_COUNT==2) && (DIM_COUNT==3)
#define TEST_BODY2DIM3
#endif

#if (BODY_COUNT==3) && (DIM_COUNT==3)
#define BODY3_DIM3
#define TEST_BODY3DIM3
#endif

#if DIM_COUNT <=0
Error : DIM_COUNT_should_be_a_positive_integer
#endif

#if BODY_COUNT <2
Error : BODY_COUNT_should_be_a_positive_integer_not_less_than_2
#endif


#endif // DEFINES_H
