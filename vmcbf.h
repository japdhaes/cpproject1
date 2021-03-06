#ifndef VMCBF_H
#define VMCBF_H
#include "vmcsolver.h"

class VMCBF:public VMCSolver
{
public:
    VMCBF();
    VMCBF(int myrank, int numprocs, int _nParticles, double _alpha, double _beta);
    void cycle(const int &i);
    void initialize();
    void setSteplength(double &sl);
    double getAcceptanceRatio(double sl);
private:
    double stepLength;
};

#endif // VMCBF_H
