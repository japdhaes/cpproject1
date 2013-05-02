#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "lib.h"
#include <armadillo>
#include <mpi.h>
#include "orbital.h"
#include "wavefunction.h"

const double pi=3.1415926535;

using namespace arma;

class VMCSolver
{
public:
    VMCSolver(int myrank, int numprocs, int _nParticles, double _alpha, double _beta);
    ~VMCSolver();

    double          distij(const mat &r, const int i, const int j);
    double          gaussianDeviate(long *seed);
    double          runMonteCarloIntegration();

    virtual void    cycle(const int &i) = 0;
    virtual void    initialize() = 0;
protected:
    int nAccepted;
    int nRejected;
    int nLocalTotalsteps;
    double deltaprobability;
    int currentparticle;
    bool closedForm;
    int numprocs, myrank;
    Wavefunction wf;
    double alpha;
    double beta;
    int nDimensions;
    int nParticles;

    double h;
    double h2;

    long idum;
    int nCycles;

    mat rOld;
    mat rNew;

    double waveFunctionOld;
    double waveFunctionNew;
    int local_nCycles;
    int argc;
    char **argv;
};

#endif // VMCSOLVER_H
