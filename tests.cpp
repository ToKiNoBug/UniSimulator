#include "tests.h"

/*
extern const double G; //gravity constant
extern const double Ms; //standard mass (solar mass)
extern const double omega_s; //standard angle speed
extern const double rho; //solar density
extern const double year; //seconds of a year
extern const double rs; //standard distance (1AU)
extern const double vs; //standard speed
extern const double as; //standard accelerate
*/

void dispConstants() {
    std::cout<<"\nG="<<G;
    std::cout<<"\nMs="<<Ms;
    std::cout<<"\nomega_s="<<omega_s;
    std::cout<<"\nrho="<<rho;
    std::cout<<"\nyear="<<year;
    std::cout<<"\nrs="<<rs;
    std::cout<<"\nvs="<<vs;
    std::cout<<"\nas="<<as<<std::endl;
}

#ifdef TEST_BODY2DIM2
void testPerformance() {
std::cerr<<"test started";

BodyVector M(Ms,1*Ms);

Statue y,dy;
y.first.setValues({{-rs,rs},
                            {0,0}});
y.second.setValues({{0,0},
                                {-vs,vs}});

DistanceMat safeMat;
Simulator::calculateSafeDistance(M,safeMat);

std::cerr<<"safe distance=\n";
std::cerr<<safeMat<<std::endl;

Interaction GM;
Simulator::calculateGM(M,GM);

clock_t curTime=std::clock();

int repeatTime=100000;

double h=1e-2*year;

for(int i=0;i<repeatTime;i++) {
    if(i%1024==0) {

        std::cerr<<i<<std::endl;

        std::cerr<<"position=\n";
        std::cerr<<y.first<<std::endl;

        std::cerr<<"velocity=\n";
        std::cerr<<y.second<<std::endl;
    }
//bool isOk=Simulator::calculateDiff(y.first,GM,safeMat,dy.second);
    bool isOk=Simulator::RK4(h,y,GM,safeMat,y);
if(!isOk) {
    std::cerr<<"stars will collide\n";
    break;
}
/*
y.first+=h*y.second;
y.second+=h*dy.second;
*/
}

std::cerr<<double(clock()-curTime)*1e6/repeatTime/CLOCKS_PER_SEC
             <<"*10^-6 s"<<std::endl;
/*
std::cerr<<"acceleration=\n";
std::cerr<<dy.second<<std::endl;
*/

//std::cerr<<(isOk?"stars won't collide\n":"stars will collide\n");

exit(0);
return;
}

void testDerivative() {
    std::cout<<"testDerivative\n";

    BodyVector M(Ms,Ms);

    Statue y;
    y.first.setValues({{-rs,rs},
                                {0,0}});
    y.second.setValues({{0,0},
                                    {-vs,vs}});

    Acceleration acc;

    DistanceMat safeMat;
    Simulator::calculateSafeDistance(M,safeMat);


    Interaction GM;
    Simulator::calculateGM(M,GM);

    std::cerr<<"\nposition:\n";
    std::cerr<<y.first;
    std::cerr<<"\nvelocity:\n";
    std::cerr<<y.second;
    Simulator::calculateDiff(y.first,GM,safeMat,acc);

    std::cerr<<"\nacceleration:\n";
    std::cerr<<acc;

    std::cout<<"\nfinished\n";
    exit(0);
    return;
}

void testSimulation() {
    std::cout<<"testSimulation\n";

    BodyVector M(Ms,Ms);

    Statue y;
    y.first.setValues({{-rs,rs},
                                {0,0}});
    y.second.setValues({{0,0},
                                    {-vs,vs}});

    Acceleration acc;

    DistanceMat safeMat;
    Simulator::calculateSafeDistance(M,safeMat);


    Interaction GM;
    Simulator::calculateGM(M,GM);

    Time maxTime=100*year;
    Time stepL=1e-2*year;

    uint64_t loopCount=maxTime/stepL;

    for(uint64_t i=0;i<loopCount;i++) {
        std::cerr<<"\nposition:\n";
        std::cerr<<y.first;
        std::cerr<<"\nvelocity:\n";
        std::cerr<<y.second;
        std::cerr<<"\nacceleration:\n";
        Simulator::calculateDiff(y.first,GM,safeMat,acc);
        std::cerr<<acc;


        bool isOk=Simulator::RK4(stepL,y,GM,safeMat,y);
        if(!isOk) {
            std::cerr<<"\nStars will collide\n";
            break;
        }
    }
    std::cout<<"\nfinished\n";
    exit(0);
    return;
}

void testEuler() {

    Simulator S;

    BodyVector mass(Ms,2*Ms);

    Statue start;
    start.first.setValues({{-rs,rs},
                                    {0,0}});
    start.second.setValues({{0,0},
                                          {-vs,vs}});

    TimeSpan tSpan=std::make_pair(0*year,10*year);

    Time step=1e-2*year;

    S.setMass(mass);

    bool noCollide;

    std::cerr<<"simulation begin\n";
    S.simulateEuler(step,tSpan,start,&noCollide);
    std::cerr<<"simulation ended\n";

    if(!noCollide) {
        std::cout<<"simulation in euler method ended by a future collision\n";
    }

    auto * sol=&S.getResult();
    for(auto it=sol->cbegin();it!=sol->cend();it++) {
        DimVector dv;
        S.calculateTotalMotion(&*it,dv);
        std::cerr<<"\nmotion:["<<dv.transpose()<<"]\n";

        double kinetic,potential;
        kinetic=S.calculateKinetic(&*it);
        potential=S.calculatePotential(&*it);

        std::cerr<<"Kinetic="<<kinetic<<" , Potential="<<potential<<std::endl;
        std::cerr<<"Energy="<<kinetic+potential<<std::endl;
    }

    exit(0);

}

void testRK4Fixed() {

    Simulator S;

    BodyVector mass(Ms,2*Ms);

    Statue start;
    start.first.setValues({{-rs,rs},
                                    {0,0}});
    start.second.setValues({{0,0},
                                          {-vs,vs}});

    TimeSpan tSpan=std::make_pair(0*year,10*year);

    Time step=0.001*year;

    S.setMass(mass);

    bool noCollide;

    std::cerr<<"simulation begin\n";
    S.simulateRK4Fixed(step,tSpan,start,&noCollide);
    std::cerr<<"simulation ended\n";

    if(!noCollide) {
        std::cout<<"simulation in euler method ended by a future collision\n";
    }

    auto * sol=&S.getResult();
    for(auto it=sol->cbegin();it!=sol->cend();it++) {
        DimVector dv;
        S.calculateTotalMotion(&*it,dv);
        std::cerr<<"\nmotion:["<<dv.transpose()<<"]\n";

        double kinetic,potential;
        kinetic=S.calculateKinetic(&*it);
        potential=S.calculatePotential(&*it);

        std::cerr<<"Kinetic="<<kinetic<<" , Potential="<<potential<<std::endl;
        std::cerr<<"Energy="<<kinetic+potential<<std::endl;
    }

    exit(0);

}

#endif
