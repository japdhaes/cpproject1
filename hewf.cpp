#include "hewf.h"

//THIS IS FOR THE HELIUM ATOM

/*It is convenient to make classes of
trial wave functions, both many-body wave functions and single-particle wave functions and the
quantum numbers involved, such as spin, orbital momentum and principal quantum numbers.*/



//this my trial wave function class
HeWF::HeWF(int myrank, int numprocs, double _alpha, double _beta):
    VMCSolver(myrank, numprocs)
{
    this->charge=2;
    this->nParticles=2;
    this->alpha=_alpha;
    this->beta=_beta;
}

double HeWF::sdratio(){
    return (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
}

double HeWF::jastrowRatio(int const k){
    return 1.0;
}

double HeWF::localEnergyClosedForm(const mat &r)
{
    double r1=norm(r.row(0),2);
    double r2=norm(r.row(1),2);
    double r12=distij(r, 0,1);
    double r1dotr2=0;
    for(int i=0; i<3; i++){
        r1dotr2+=r(0,i)*r(1,i);
    }


    double EL1 = (this->alpha-this->charge)*(1.0/r1+1/r2)+1.0/r12-this->alpha*this->alpha;
    double braces=(this->alpha*(r1+r2)/r12*(1-r1dotr2/r1/r2)-1.0/2/(1+this->beta*r12)/(1+this->beta*r12) - 2.0/r12 + 2*this->beta/(1+this->beta*r12));
    double EL2=EL1+1.0/2/(1+this->beta*r12)/(1+this->beta*r12)*braces;

    return EL2;
}

double HeWF::waveFunction(const mat &r)
{
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
        }
    }
    r12=sqrt(r12);

    double argument = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }
    return exp(-argument * alpha)*exp(r12/2/(1+beta*r12));
}
