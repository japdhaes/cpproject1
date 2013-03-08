#ifndef BEWF_H
#define BEWF_H

#include "vmcsolver.h"

class BeWF : public VMCSolver
{
public:
    BeWF();
    BeWF(int myrank, int numprocs, double _alpha, double _beta);

    double localEnergyClosedForm(const mat &r);
    double waveFunction(const mat &r);
    double aSDWF(const mat &r);
    double wavefunction(const mat &r, const mat &fij);
    double jastrowWF(const mat &fij);
    double totalSD(const mat &r);

    double jastrowRatio(const int k);
    double sdratio();

//    double (BeWF::*hydrogenWF[2]) (const vec3 &position);
};

#endif // BEWF_H
