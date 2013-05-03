#include "vmcbf.h"


VMCBF::VMCBF(int myrank, int numprocs, int _nParticles, double _alpha, double _beta):
    VMCSolver(myrank,numprocs,_nParticles, _alpha, _beta),
    stepLength(1.0)
{

}

void VMCBF::cycle(const int &i){
    for(int j = 0; j < nDimensions; j++) {
        rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
    }
    this->wf.setNewPos(rNew);

    // Recalculate the value of the wave function
    waveFunctionNew = wf.evaluate(rNew);
//    cout << waveFunctionNew <<endl;

    //accepting
    if(ran2(&idum) <= wf.calcRatio()) {
        wf.acceptMove();
        rOld.row(i)=rNew.row(i);
    }
    //rejecting
    else {
        wf.rejectMove();
        rNew.row(i)=rOld.row(i);
    }
}

void VMCBF::initialize()
{
    rNew = rOld;
    wf.initialize(rOld);
}
