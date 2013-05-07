#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "lib.h"
#include <armadillo>
#include <mpi.h>
#include "orbital.h"
#include "wavefunction.h"
#include "datalogger.h"

const double pi=3.1415926535;
const int nrOfCyclesEachOutput=1e4;

using namespace arma;

class VMCSolver
{
public:
    VMCSolver(int myrank, int numprocs, int _nParticles, double _alpha, double _beta);
    ~VMCSolver();

    double          gaussianDeviate(long *seed);
    double          runMonteCarloIntegration();

    virtual void    cycle(const int &i) = 0;
    virtual void    initialize() = 0;
    void            solverInitializer();

    void            setCycles(const int &_nCycles);


    void            setOutput(bool output);
    void            setOutputDirectory(const string &directoryName);

protected:
    int             nAccepted;
    int             nRejected;
    int             nDimensions;
    int             nParticles;

    double          deltaprobability;
    int             currentparticle;
    bool            closedForm;

    double          h;
    double          h2;

    long            idum;
    Wavefunction    wf;

    mat             rOld;
    mat             rNew;

    double          waveFunctionOld;
    double          waveFunctionNew;

private:
    int             local_nCycles;
    int             nLocalTotalsteps;
    int             numprocs, myrank;
    int             nCycles;


    bool            createOutput;
    string          outputDirectory;
    string          outputFilename;
    Datalogger      datalogger;

};

#endif // VMCSOLVER_H
