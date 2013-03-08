#ifndef WAFEFUNCTION_H
#define WAFEFUNCTION_H

#include "vmcsolver.h"


class HeWF : public VMCSolver
{
public:
    HeWF();
    HeWF(int myrank, int numprocs, double _alpha, double _beta);

    double jastrowRatio(const int k);
    double sdratio();

    double localEnergyClosedForm(const state &astate);
    double waveFunction(const state &astate);
};

#endif // WAFEFUNCTION_H
