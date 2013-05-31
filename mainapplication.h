#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

#include <iostream>
#include "lib.h"
#include <ctime>
#include "vmcbf.h"
#include "vmcis.h"
#include <unistd.h>
#include "testdirectory.h"
#include "minimizer.h"

class MainApplication
{
public:
    MainApplication(int _myrank, int _numprocs);
    void runApplication();
    void runSimulation();
    void simulateWithOutput(int nParticles, double alpha, double beta);
    void runSimulation(bool importancesampling, int nCycles, int nParticles, double alpha, double beta);
    void minimizeBruteForce();
    void minimizeNM();
    void steplengthSecant();
private:
    int numprocs, myrank;
    double alphaHe;
    double betaHe;
    double alphaBe;
    double betaBe;
    double alphaNe;
    double betaNe;
    double alphaH2;
    double betaH2;
    double distH2;
};

#endif // MAINAPPLICATION_H
