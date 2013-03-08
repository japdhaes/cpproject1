#include "bewf.h"

//THIS IS FOR THE BERYLLIUM ATOM

/*It is convenient to make classes of
trial wave functions, both many-body wave functions and single-particle wave functions and the
quantum numbers involved, such as spin, orbital momentum and principal quantum numbers.*/


//this my trial wave function class
BeWF::BeWF(int myrank, int numprocs, double _alpha, double _beta):
    VMCSolver(myrank, numprocs)
{
    this->charge=4;
    this->nParticles=4;
    this->alpha=_alpha;
    this->beta=_beta;

//    this->hydrogenWF[0]=&BeWF::phi1s2;
//    this->hydrogenWF[1]=&BeWF::phi2s2;
}

double BeWF::localEnergyClosedForm(const state &astate)
{
    cout << __LINE__<<endl;
    //impossible
    return 0;
}

double BeWF::waveFunction(const state &astate){
    const mat r = astate.r;
//    mat fij = calculatefij()
      return totalSD(astate);
//    return this->aSDWF(r);
}


double BeWF::wavefunction(const state &astate){
    return aSDWF(astate)*jastrowWF(astate);
}

double BeWF::jastrowWF(const state &astate){
    const mat fij=astate.fij;
    double arg=0.0;
    for(int i=0; i<nParticles; i++){
        for(int j=i+1; j<nParticles;j++){

            arg+=fij(i,j);
        }
    }
    return exp(arg);
}

double BeWF::totalSD(const state &astate){
    const mat r = astate.r;
    int halfparticles=nParticles/2;
    mat sdup = zeros<mat>(halfparticles, halfparticles);
    mat sddown = zeros<mat>(halfparticles, halfparticles);
    vec3 temp;
    for(int i=0; i<halfparticles; i++){
        for(int j=0; j<halfparticles;j++){
            temp = r.row(j);
            sdup(i,j)=this->hydrogenWF(i, temp);
            temp = r.row(j+halfparticles);
            sddown(i,j)=this->hydrogenWF(i, temp);
        }
    }
    return det(sdup)*det(sddown);
}

double BeWF::sdratio(){
    return (this->newS.wavefunction*this->newS.wavefunction) / (this->oldS.wavefunction*this->oldS.wavefunction);
}

double BeWF::jastrowRatio(int const k){

    return 1.0;
    double dU=0.0;
    for(int i=0; i<k;i++){
        dU+=calcElementfij(this->newS,i,k)-calcElementfij(this->oldS,i,k);
    }
    for(int i=k+1; i<nParticles;i++){
        dU+=calcElementfij(this->newS,k,i)-calcElementfij(this->oldS,k,i);
    }
    return exp(dU);
}



//explicitly the slater determinant split up between spin up and spin down particles
//explicitly calculating determinant of the split up parts
double BeWF::aSDWF(const state &astate)
{
//    cout << r<< endl;
    const mat r=astate.r;
    double dist[nParticles];
    for(int i=0; i<nParticles;i++){
        dist[i]=0.0;
        for(int j=0; j<nDimensions; j++){
            dist[i]+=r(i,j)*r(i,j);
        }
        dist[i]=sqrt(dist[i]);
    }
    double answer=1.0;
    for(int i=0; i<4; i+=2){
        answer*=(phi1s(dist[i])*phi2s(dist[i+1])-phi1s(dist[i+1])*phi2s(dist[i]));
    }
    return answer;
}



