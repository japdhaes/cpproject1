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

    int local_nCycles;
    int argc;
    char **argv;
    double distij(const mat &r, const int i, const int j);
    int nAccepted;
    int nRejected;
    int nLocalTotalsteps;
    double deltaprobability;

    bool importanceSampling;
    bool closedForm;

    double D;
    mat quantumForce(const mat &r, const double &wavefunction);
    int numprocs, myrank;

    double hydrogenWF(const int &i, vec3 &r);
    double phi1s(const vec3 &position);
    double phi2s(const vec3 &position);
    double phi1s(const double position);
    double phi2s(const double position);



    virtual double localEnergyClosedForm(const mat &r)=0;
    virtual double waveFunction(const mat &r)=0;

    void calculaterij(const mat &r, mat &rij);
    void updaterij(const mat &r, mat &rij, const int j);
    virtual double jastrowRatio(const int k)=0;
    virtual double sdratio()=0;

    void cycleIS(const int &i);
    void cycleWithoutIS(const int &i);
    double gaussianDeviate(long *seed);
    double calcElementfij(const mat &rij, const int i, const int j);
    void calculatefij(const mat &rij, mat &fij);
protected:
    double alpha;
    double beta;

    double localEnergy(const mat &r);

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
