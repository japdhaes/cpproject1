#ifndef BEWF_H
#define BEWF_H

#include "vmcsolver.h"

class BeWF : public VMCSolver
{
public:
    BeWF();
    BeWF(int myrank, int numprocs, double _alpha, double _beta);

    double localEnergyClosedForm(const state &astate);
    double waveFunction(const state &astate);
    double aSDWF(const state &astate);
    double wavefunction(const state &astate);
    double jastrowWF(const state &astate);
    double totalSD(const state &astate);

    double jastrowRatio(const int k);
    double sdratio();

    //    double (BeWF::*hydrogenWF[2]) (const vec3 &position);

};

#endif // BEWF_H
