#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "lib.h"
#include <armadillo>
#include <mpi.h>

const double pi=3.1415926535;

using namespace arma;

class VMCSolver
{
public:
    VMCSolver(int myrank, int numprocs);
    ~VMCSolver();
    double runMonteCarloIntegration();
    struct state{
        double wavefunction;
        mat r, rij, qForce, fij;
    };

    int local_nCycles;
    int argc;
    char **argv;
    double distij(const state &astate, const int &i, const int &j);
    int nAccepted;
    int nRejected;
    int nLocalTotalsteps;
    double deltaprobability;

    bool importanceSampling;
    bool closedForm;

    double D;
    mat quantumForce(const state &astate);
    int numprocs, myrank;

    double hydrogenWF(const int &i, vec3 &r);
    double phi1s(const vec3 &position);
    double phi2s(const vec3 &position);
    double phi1s(const double position);
    double phi2s(const double position);



    virtual double localEnergyClosedForm(const state &astate)=0;
    virtual double waveFunction(const state &astate)=0;

    void calculaterij(state &astate);
    void updaterij(state &astate, const int &j);
    virtual double jastrowRatio(const int k)=0;
    virtual double sdratio()=0;

    void cycleIS(const int &i);
    void cycleWithoutIS(const int &i);
    double gaussianDeviate(long *seed);
    void calculatefij(state &astate);


    state newS;
    state oldS;
    double calcElementfij(const state &astate, const int &i, const int &j);


protected:
    double alpha;
    double beta;

    double localEnergy(const state &astate);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    int nCycles;

    mat rOld;
    mat rNew;

    mat rijOld;
    mat rijNew;

    mat qForceOld;
    mat qForceNew;

    double waveFunctionOld;
    double waveFunctionNew;

    mat fijOld;
    mat fijNew;
};

#endif // VMCSOLVER_H
