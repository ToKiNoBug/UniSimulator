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

#ifdef IS_POINT_TUPLE
const uint8_t tupleTimeIndex=0,
                                tuplePositionIndex=1,
                                tupleVelocityIndex=2,
                                tupleAccelerlationIndex=3;
#endif

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
if(!isOk) return false;

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

    for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
        Eigen::Tensor<double,0> temp=y.first.chip(dim,0).mean();
        double && meanPos=std::move(temp(0));

        temp=y.second.chip(dim,0).mean();
        double && meanVelocity=std::move(temp(0));

        for(uint16_t body=0;body<BODY_COUNT;body++) {
            y.first(dim,body)-=meanPos;
            y.second(dim,body-=meanVelocity);
        }

    }

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

    for(uint16_t dim=0;dim<DIM_COUNT;dim++) {
        Eigen::Tensor<double,0> temp=y.first.chip(dim,0).mean();
        double && meanPos=std::move(temp(0));

        temp=y.second.chip(dim,0).mean();
        double && meanVelocity=std::move(temp(0));

        for(uint16_t body=0;body<BODY_COUNT;body++) {
            y.first(dim,body)-=meanPos;
            y.second(dim,body-=meanVelocity);
        }

    }

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


void Simulator::setMass(const MassVector & _mass) {
    mass=_mass;
}

const MassVector & Simulator::getMass() const {
    return mass;
}

const std::list<Point> & Simulator::getResult() const {
    return sol;
}

double Simulator::calculateKinetic(const std::_List_const_iterator<Point> it) const {
    Eigen::Tensor<double,1> speedSquare
            =it->second.second.square().sum(Eigen::array<int,1>({0}));
    MassVector v;
    for(uint16_t i=0;i<BODY_COUNT;i++) {
        v(i)=speedSquare(i)/2;
    }

    return (v*mass).sum();
}

double Simulator::calculatePotential(const std::_List_const_iterator<Point> it) const {
    DistanceMat realDistance;

    realDistance.setZero();

    const auto & pos=it->second.first;

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

    return realDistance.sum();
}

double Simulator::calculateEnergy(const std::_List_const_iterator<Point> it) const {
    return calculateKinetic(it)+calculatePotential(it);
}

void Simulator::calculateTotalMotion(const std::_List_const_iterator<Point> it,
                          DimVector & dest) const {
Eigen::Array<double,DIM_COUNT,BODY_COUNT> speed;
for(uint16_t dim=0;dim<DIM_COUNT;dim++)
    for(uint16_t body=0;body<BODY_COUNT;body++) {
        speed(dim,body)=it->second.second(dim,body);
    }

speed.rowwise()*=mass.transpose();

dest=speed.rowwise().sum();
}

#if DIM_COUNT <=0
DIM_COUNT_should_be_a_positive_integer
#endif

#if BODY_COUNT <2
BODY_COUNT_should_be_a_positive_integer_no_less_than_2
#endif
