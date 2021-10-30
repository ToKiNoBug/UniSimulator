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

#include "Simulator.h"

//const uint32_t BodyCount=3;
//const Eigen::Index DimCount=3;
const double G=6.67259e-11;
const double Ms=2e30;
const double year=365*24*60*60;
const double omega_s=2*M_PI/year;
const double rho=1.409e3;
const double rs=0.5*std::pow(2*G*Ms/(omega_s*omega_s),1.0/3);
const double vs=omega_s*rs;
const double as=omega_s*vs;
const double INF=1e300;

#ifdef IS_POINT_TUPLE
const uint8_t tupleTimeIndex=0,
                                tuplePositionIndex=1,
                                tupleVelocityIndex=2,
                                tupleAccelerlationIndex=3;
#endif

Simulator::Simulator() {
sol.clear();
}

void Simulator::calculateSafeDistance(const BodyVector & mass,
                                                                DistanceMat &dest) {
    //dest.setZero();

    BodyVector radius=(3*mass/(4*M_PI*rho)).pow(1.0/3);
    auto massMat=radius.replicate(1,BODY_COUNT);

    //dest=2.44*massMat.max(massMat.transpose());
    dest=2.44*massMat.max(massMat.transpose());
}

void Simulator::calculateGM(const BodyVector & mass,
                           Interaction & GM) {

for(int i=0;i<GM.dimension(1);i++) {

    GM.chip(i,1).setConstant(G*mass(i));

}
//std::cerr<<"caculate finished\n";
}

bool Simulator::calculateDiff(const Position & pos,
                             const Interaction & GM,
                             const DistanceMat & safeDistance,
                             Acceleration & dy) {
Interaction T;
T.setZero();

DistanceMat realDistance;

realDistance.setConstant(INF);

//std::cerr<<"\npos.chip(0,1)=\n"<<pos.chip(0,1);

for(uint16_t i=0;i<BODY_COUNT;i++) {
    for(uint16_t j=i+1;j<BODY_COUNT;j++) {

        auto delta_r=pos.chip(i,1)-pos.chip(j,1);
        Eigen::Tensor<double,0> distanceT=delta_r.square().sum();
        double && distance= std::move(distanceT(0));
        double && distanceSqrt=std::sqrt(distance);

        realDistance(i,j)=distanceSqrt;
        realDistance(j,i)=distanceSqrt;

        T.chip(j,2).chip(i,1)=delta_r/(distanceSqrt*distance);
        T.chip(i,2).chip(j,1)=-T.chip(j,2).chip(i,1);
    }
}

//std::cerr<<"\nrealDistance:\n"<<realDistance;

if((realDistance<safeDistance).any()) {
    return false;
}

T*=GM;
//std::cerr<<"\nT=\n"<<T;
//dy.first=y.second;

dy=T.sum(Eigen::array<int,1>({1}));
return true;
}

bool Simulator::RK4(const Time h,
                    const Statue & y_n,
                    const Interaction & GM,
                    const DistanceMat & safeDistance,
         Statue & y_n1) {
bool isOk=true;

Time && halfStep=h/2;

Derivative k1,k2,k3,k4;
Position temp;

k1.first=y_n.second;
isOk=calculateDiff(y_n.first,GM,safeDistance,k1.second);
if(!isOk) {
    //y_n1.second.setConstant(INF);
    //make error extremely obvious when collide happens
    return false;
}

k2.first=y_n.second+halfStep*k1.second;
temp=y_n.first+halfStep*k1.first;
isOk=calculateDiff(temp,GM,safeDistance,k2.second);

k3.first=y_n.second+halfStep*k2.second;
temp=y_n.first+halfStep*k2.first;
isOk=calculateDiff(temp,GM,safeDistance,k3.second);

k4.first=y_n.second+h*k3.second;
temp=y_n.first+h*k3.first;
isOk=calculateDiff(temp,GM,safeDistance,k4.second);

y_n1.first=y_n.first+h/6*(k1.first+2*k2.first+2*k3.first+k4.first);
y_n1.second=y_n.second+h/6*(k1.second+2*k2.second+2*k3.second+k4.second);

return true;
}

void Simulator::clear() {
    sol.clear();
}

void Simulator::simulateEuler(const double step,
                   TimeSpan tSpan,
                   Statue y,
                   bool * noCollide) {
    if(tSpan.second==tSpan.first) {
        std::cerr<<"the begging time shouldn't equal to end time\n";
        return;
    }

    if(tSpan.second<tSpan.first) {
        std::swap(tSpan.first,tSpan.second);
    }

    if(tSpan.second-tSpan.first<step) {
        std::cerr<<"time span should be greater than time step\n";
        return;
    }

    positonAlign(y.first);
    motionAlign(mass,y.second);

    clear();

    DistanceMat safeDistance;
    calculateSafeDistance(mass,safeDistance);

    Interaction GM;
    calculateGM(mass,GM);

    Acceleration acc;

    uint64_t i=0;
    Time curTime;

    if(noCollide!=nullptr) {
        *noCollide=true;
    }

    while (true) {
        curTime=tSpan.first+i*step;

        if(curTime>tSpan.second) {
            break;
        }

        sol.push_back(std::make_pair(curTime,y));

        bool isOk=calculateDiff(y.first,GM,safeDistance,acc);

        if(!isOk) {
            std::cerr<<"Stars will collide\n";
            if(noCollide!=nullptr) {
                *noCollide=false;
            }
            break;
        }

        y.first+=step*y.second;
        y.second+=step*acc;
        i++;
    }
}

void Simulator::simulateRK4Fixed(const double step,
                   TimeSpan tSpan,
                   Statue y,
                   bool * noCollide) {
    if(tSpan.second==tSpan.first) {
        std::cerr<<"the begging time shouldn't equal to end time\n";
        return;
    }

    if(tSpan.second<tSpan.first) {
        std::swap(tSpan.first,tSpan.second);
    }

    if(tSpan.second-tSpan.first<step) {
        std::cerr<<"time span should be greater than time step\n";
        return;
    }

    clear();


    positonAlign(y.first);
    motionAlign(mass,y.second);

    DistanceMat safeDistance;
    calculateSafeDistance(mass,safeDistance);

    Interaction GM;
    calculateGM(mass,GM);

    Acceleration acc;

    uint64_t i=0;
    Time curTime;

    if(noCollide!=nullptr) {
        *noCollide=true;
    }

    while (true) {
        curTime=tSpan.first+i*step;

        if(curTime>tSpan.second+step) {
            break;
        }

        sol.push_back(std::make_pair(curTime,y));

        bool isOk=RK4(step,y,GM,safeDistance,y);

        if(!isOk) {
            std::cerr<<"Stars will collide\n";
            if(noCollide!=nullptr) {
                *noCollide=false;
            }
            break;
        }
        i++;
    }
}

void Simulator::simulateRK4Var1(double step,
                                TimeSpan tSpan, Statue y, bool *noCollide) {
    if(tSpan.second==tSpan.first) {
        std::cerr<<"the begging time shouldn't equal to end time\n";
        return;
    }

    if(tSpan.second<tSpan.first) {
        std::swap(tSpan.first,tSpan.second);
    }

    if(tSpan.second-tSpan.first<step) {
        std::cerr<<"time span should be greater than time step\n";
        return;
    }

    clear();

    positonAlign(y.first);
    motionAlign(mass,y.second);

    DistanceMat safeDistance;
    calculateSafeDistance(mass,safeDistance);

    Interaction GM;
    calculateGM(mass,GM);

    Acceleration acc;
    Time curTime=tSpan.first;

    Statue y_h,y_h_2;

    static const double searchRatio=0.8;
    static const int rank=4;
    static const double ratio=std::pow(2,rank)-1;

    while(true) {
        sol.push_back(std::make_pair(curTime,y));

        if(curTime>tSpan.second) {
            break;
        }

        double minStep=16*std::nextafter(curTime,curTime*2)-16*curTime;

        bool isOk=true;

        isOk=RK4(step,y,GM,safeDistance,y_h);
        if(!isOk) {
            if(noCollide!=nullptr) {
                *noCollide=isOk;
            }
            std::cerr<<"stars will collide (Line"<<__LINE__<<")\n";
            break;
        }

        RK4(step/2,y,GM,safeDistance,y_h_2);

        RK4(step/2,y_h_2,GM,safeDistance,y_h_2);

        if(isErrorTolerantable(y_h,y_h_2)) {
            //error is tolerantable, scale up until next value is untolerantable
            while(true) {
                step/=searchRatio;
                //std::cerr<<"step increasing to "<<step<<std::endl;
                RK4(step,y,GM,safeDistance,y_h);
                RK4(step/2,y,GM,safeDistance,y_h_2);
                RK4(step/2,y_h_2,GM,safeDistance,y_h_2);
                if(!isErrorTolerantable(y_h,y_h_2)) {
                    step*=searchRatio;
                    RK4(step,y,GM,safeDistance,y_h);
                    RK4(step/2,y,GM,safeDistance,y_h_2);
                    RK4(step/2,y_h_2,GM,safeDistance,y_h_2);
                    break;
                }
            }
        } else {
            while(true) {
                step*=searchRatio;
                //std::cerr<<"step shrinking to "<<step<<std::endl;
                if(step<=minStep) {
                    if(noCollide!=nullptr) {
                        *noCollide=false;
                    }
                    std::cerr<<"stars will collide\n";
                    break;
                }
                RK4(step,y,GM,safeDistance,y_h);
                RK4(step/2,y,GM,safeDistance,y_h_2);
                RK4(step/2,y_h_2,GM,safeDistance,y_h_2);
                if(isErrorTolerantable(y_h,y_h_2)) {
                    break;
                }
            }
        }
        //std::cerr<<"step accepted ="<<step<<std::endl;
        Statue yh2_error;
        yh2_error.first=(y_h_2.first-y_h.first)/ratio;
        yh2_error.second=(y_h_2.second-y_h.second)/ratio;

        y.first=y_h_2.first+yh2_error.first;
        y.second=y_h_2.second+yh2_error.second;
        curTime+=step;

    }

}

bool Simulator::isErrorTolerantable(const Statue &y_h,
                                    const Statue &y_h_2,
                                    double errorRatio) {
    static const int rank=4;
    static const double ratio=std::pow(2,rank)-1;

    errorRatio*=ratio;

    bool isTolerantable=true;
    auto posError=(y_h_2.first-y_h.first).abs();

    Eigen::Tensor<bool,0> temp=(posError>=(errorRatio*rs)).any();
    isTolerantable=!temp(0);
    if(!isTolerantable)
        return false;
    auto velocityError=(y_h_2.second-y_h.second).abs();
    temp=(velocityError>=(errorRatio*vs)).any();
    isTolerantable=!temp(0);

    return isTolerantable;
}

void Simulator::setMass(const BodyVector & _mass) {
    mass=_mass;
}

const BodyVector & Simulator::getMass() const {
    return mass;
}

const std::list<Point> & Simulator::getResult() const {
    return sol;
}

double Simulator::calculateKinetic(const Statue & it) const {
    Eigen::Tensor<double,1> speedSquare
            =it.second.square().sum(Eigen::array<int,1>({0}));
    BodyVector v;
    for(uint16_t i=0;i<BODY_COUNT;i++) {
        v(i)=speedSquare(i)/2;
    }

    return (v*mass).sum();
}

double Simulator::calculatePotential(const Statue & it) const {
    DistanceMat realDistance;

    realDistance.setZero();

    const auto & pos=it.first;

    for(uint16_t i=0;i<BODY_COUNT;i++) {
        for(uint16_t j=i+1;j<BODY_COUNT;j++) {

            auto delta_r=pos.chip(i,1)-pos.chip(j,1);
            Eigen::Tensor<double,0> distanceT=delta_r.square().sum();
            //double && distance= std::move(distanceT(0));
            double && distanceSqrt=-G/std::sqrt(distanceT(0));

            realDistance(i,j)=distanceSqrt;
            realDistance(j,i)=distanceSqrt;
        }
    }
    auto massMatrix=mass.matrix();
    realDistance*=(massMatrix*massMatrix.transpose()).array();

    return realDistance.sum()/2;
}

double Simulator::calculateEnergy(const Statue & it) const {
    return calculateKinetic(it)+calculatePotential(it);
}

void Simulator::calculateTotalMotion(const Statue & it,
                          DimVector & dest) const {
Eigen::Array<double,DIM_COUNT,BODY_COUNT> speed;
for(uint16_t dim=0;dim<DIM_COUNT;dim++)
    for(uint16_t body=0;body<BODY_COUNT;body++) {
        speed(dim,body)=it.second(dim,body);
    }

speed.rowwise()*=mass.transpose();

dest=speed.rowwise().sum();
}

void Simulator::deval(const Simulator *source,
                          Simulator *dest,
                          const Eigen::ArrayXd & timeQueried) {
if(source->sol.back().first<timeQueried.maxCoeff()) {
    std::cerr<<"Error when interploting! outerplot is invalid\n";
    std::cerr<<"source ended at time "<<source->sol.back().first
            <<" but quried at time "<<timeQueried.maxCoeff()<<std::endl;
    return;
}

if(source->sol.size()<=1) {
    std::cerr<<"source's sol has fewer than 2 statues, can't interplot in such condition\n";
    return;
}

dest->clear();
dest->mass=source->mass;

Interaction GM;
DistanceMat safeDistance;
calculateGM(source->mass,GM);
calculateSafeDistance(source->mass,safeDistance);

//assume that timeQueried is sorted into ascending order
auto cur=source->sol.cbegin();
auto next=cur;      next++;


Point curP=*cur;
    for(auto curTime : timeQueried)
    {
        while(curTime>next->first) {
            cur++;
            next++;
            curP=*cur;
            //std::cout<<"region marched a step\n";
        }
        /*
        if(next==source->sol.cend()) {
            break;
        }*/
        //std::cout<<"interplot time "<<curTime<<" in range ["<<cur->first<<" , "<<next->first<<"]\n";
        double step=curTime-curP.first;
        Simulator::RK4(step,curP.second,GM,safeDistance,curP.second);
        curP.first=curTime;
        dest->sol.push_back(curP);
    }

}

void Simulator::motionAlign(const BodyVector &mass, Velocity &velocity) {
    Eigen::Tensor<double,1> massTensor(mass.size());

    for(uint32_t body=0;body<BODY_COUNT;body++) {
        massTensor[body]=mass[body];
    }

    double massSum=mass.sum();

    for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
        Eigen::Tensor<double,0> temp=(velocity.chip(dim,0)*massTensor).sum();
        double && deltaVelocity=temp(0)/massSum;

        for(uint16_t body=0;body<BODY_COUNT;body++) {
            velocity(dim,body)-=deltaVelocity;
        }

    }
}

void Simulator::positonAlign(Position & pos) {
    for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
        Eigen::Tensor<double,0> temp=pos.chip(dim,0).mean();
        double && meanPos=std::move(temp(0));

        for(uint32_t body=0;body<BODY_COUNT;body++) {
            pos(dim,body)-=meanPos;
        }
    }
}
