#include "vmcis.h"
#include "zignor.c"
#include "zigrandom.c"

VMCIS::VMCIS(int myrank, int numprocs, int _nParticles, double _alpha, double _beta):
    VMCSolver(myrank,numprocs,nParticles, alpha, beta),
    D(0.5)
{

}

void VMCIS::cycle(const int &i){
    for(int j = 0; j < nDimensions; j++) {
        rNew(i,j) = rOld(i,j) + this->gaussianDeviate(&idum)*sqrt(h)+qForceOld(i,j)*D*h;
    }
    this->wf.setNewPos(rNew);

    // Recalculate the value of the wave function
    waveFunctionNew = wf.evaluate(rNew);
    qForceNew=this->quantumForce(rNew, waveFunctionNew);

    // Check for step acceptance (if yes, update position, if no, reset position)
    double expratio = 0.0;
    for(int j=0;j<nDimensions;j++){
        expratio+=(qForceOld(i,j)+qForceNew(i,j))*(D*h*(qForceOld(i,j)-qForceNew(i,j))+2.0*(rOld(i,j)-rNew(i,j)));
    }
    expratio/=4.0;
    expratio=exp(expratio);

    //accepting
    if(ran2(&idum) <= expratio*wf.getRatio()) {
        wf.acceptMove();
        rOld.row(i)=rNew.row(i);
        nAccepted++;
    }
    //rejecting
    else {
        wf.rejectMove();
        rNew.row(i)=rOld.row(i);
        nRejected++;
    }
}

mat VMCIS::quantumForce(const mat &r, const double &wavefunction){
    mat qforce(nParticles, nDimensions);
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);
    double waveFunctionMinus, waveFunctionPlus;
    rPlus = rMinus = r;

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = wf.evaluate(rMinus);
            waveFunctionPlus = wf.evaluate(rPlus);
            qforce(i,j) = waveFunctionPlus-waveFunctionMinus;
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    qforce/=(wavefunction*h);
    return qforce;
}

void VMCIS::initialize()
{
    rNew = rOld;
    wf.initialize(rOld);
    waveFunctionOld = wf.evaluate(rOld);
    qForceNew = qForceOld = quantumForce(rOld, waveFunctionOld) ;
}
