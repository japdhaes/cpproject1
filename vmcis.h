#ifndef VMCIS_H
#define VMCIS_H
#include "vmcsolver.h"

class VMCIS:public VMCSolver
{
public:
    VMCIS(int myrank, int numprocs, int _nParticles, double _alpha, double _beta);
    mat quantumForce();
    void initialize();
    void cycle(const int &i);
private:
    double D;
    mat qForceOld;
    mat qForceNew;
};

#endif // VMCIS_H
