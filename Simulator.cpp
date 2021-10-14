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
const double INF=1e300;

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

void Simulator::calculateGM(const MassVector & mass,
                           Interaction & GM) {


    //GM.setRandom();

    //std::cerr<<"the 2nd dim has "<<GM.dimension(1)<<" dims\n";
    //auto dims=GM.dimensions();
    //std::cerr<<"size of GM is "<<dims[0]<<" , "<<dims[1]<<" , "<<dims[2]<<"\n";


    //int j=1;
    //Eigen::Tensor<double,2> subGM=GM.chip(j,1);

    //std::cerr<<"GM=\n";
    //std::cerr<<GM<<std::endl;

    //std::cerr<<"subGM=\n";
    //std::cerr<<subGM<<std::endl;

for(int i=0;i<GM.dimension(1);i++) {
    //std::cerr<<"i="<<i<<std::endl;
    //std::cerr<<GM.chip(i,1)<<std::endl;
    GM.chip(i,1).setConstant(G*mass(i));
    //std::cerr<<GM.chip(i,1)<<std::endl;
}
//std::cerr<<"caculate finished\n";
}

bool Simulator::calculateDiff(const Statue & y,
                             const Interaction & GM,
                             const DistanceMat & safeDistance,
                             Statue & dy) {
Interaction T;
T.setZero();

DistanceMat realDistance;

realDistance.setConstant(INF);

for(uint16_t i=0;i<BODY_COUNT;i++) {
    for(uint16_t j=i+1;j<BODY_COUNT;j++) {

        auto delta_r=y.first.chip(i,1)-y.first.chip(j,1);
        Eigen::Tensor<double,0> distanceT=delta_r.square().sum();
        double && distance= std::move(distanceT(0));

        realDistance(i,j)=distance;
        realDistance(j,i)=distance;

        T.chip(j,2).chip(i,1)=delta_r/(std::sqrt(distance)*distance);
        T.chip(i,2).chip(j,1)=-T.chip(j,2).chip(i,1);
    }
}

if((realDistance<safeDistance).any()) {
    return false;
}

T*=GM;

//dy.first=y.second;

dy.second=T.sum(Eigen::array<int,1>({1}));
return true;
}


#if DIM_COUNT <=0
DIM_COUNT_should_be_a_positive_integer
#endif

#if BODY_COUNT <2
BODY_COUNT_should_be_a_positive_integer_no_less_than_2
#endif
