#include "Simulator.h"

//const uint32_t BodyCount=3;
//const Eigen::Index DimCount=3;
const double G=6.67259e-11;
const double Ms=2e30;
const double omega_s=2e-7;
const double rho=1.409e3;
const double year=2*M_PI/omega_s;
const double rs=0.5*std::pow(2*G*Ms/(omega_s*omega_s),1.0/3);
const double vs=omega_s*rs;
const double as=omega_s*vs;
const double TimeMax=300*year;

Simulator::Simulator() {
sol.clear();
}

void Simulator::calculateSafeDistance(const MassVector & mass,
                                                                DistanceMat &dest) {
    //dest.setZero();

    MassVector radius=(3*mass/(4*M_PI*rho)).pow(1.0/3);
    auto massMat=radius.replicate(1,BODY_COUNT);

    dest=2.44*massMat.max(massMat.transpose());
}

bool Simulator::calculateDiff(const Statue & y,
                             const Interaction & GM,
                             const DistanceMat & safeDistance,
                             Statue & dy) {
Interaction T;
T.setZero();

DistanceMat realDistance;

realDistance.setZero();

for(uint16_t i=0;i<BODY_COUNT;i++) {
    for(uint16_t j=i+1;j<BODY_COUNT;j++) {
        Eigen::TensorFixedSize<double,Eigen::Sizes<DIM_COUNT>>
                delta_r=y.first.chip(i,0)-y.first.chip(j,0);
        Eigen::TensorFixedSize<double,Eigen::Sizes<1>> distanceT=delta_r.square().sum();
        double && distance= std::move(distanceT(0));

        realDistance(i,j)=distance;
        realDistance(j,i)=distance;

        T.chip(j,2).chip(i,1)=delta_r/(std::sqrt(distance)*distance);
    }
}

if((realDistance<safeDistance).any()) {
    return false;
}

T*=GM;

//dy.first=y.second;

auto temp=T.sum(Eigen::array<int,1>({1}));

std::cerr<<"size of T.sum(Eigen::array<int,1>({1}))=["<<temp.NumDimensions<<std::endl;

dy.second=temp;
return true;
}
