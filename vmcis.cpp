#include "vmcis.h"

VMCIS::VMCIS(int myrank, int numprocs, int _nParticles, double _alpha, double _beta):
    VMCSolver(myrank,numprocs,_nParticles, _alpha, _beta),
    D(0.5)
{

}


void VMCIS::cycle(const int &i){
    for(int j = 0; j < nDimensions; j++) {
        rNew(i,j) = rOld(i,j) + this->gaussianDeviate(&idum)*sqrt(h)+qForceOld(i,j)*D*h;
    }
    this->wf.setNewPos(rNew);

    // Calculate new quantum force
    qForceNew=this->quantumForce();

    // Check for step acceptance (if yes, update position, if no, reset position)
    double greensratio = 0.0;
    for(int k=0; k<nParticles;k++){
        for(int j=0;j<nDimensions;j++){
            greensratio+=(qForceOld(k,j)+qForceNew(k,j))*(D*h*(qForceOld(k,j)-qForceNew(k,j))+2.0*(rOld(k,j)-rNew(k,j)));
        }
    }
    greensratio/=4.0;
    greensratio=exp(greensratio);

    //accepting
    if(ran2(&idum) <= greensratio*wf.calcRatio()) {
        wf.acceptMove();
        rOld.row(i)=rNew.row(i);
        qForceOld=qForceNew;
        nAccepted++;
    }
    //rejecting
    else {
        wf.rejectMove();
        rNew.row(i)=rOld.row(i);
        qForceNew=qForceOld;
        nRejected++;
    }
}

mat VMCIS::quantumForce(){
    mat qforce(nParticles, nDimensions);
    qforce=2.0*wf.localGradient();
    return qforce;
}

void VMCIS::initialize()
{
    solverInitializer();
    rNew = rOld;
    wf.initialize(rOld);
    qForceNew = qForceOld = quantumForce() ;
}
