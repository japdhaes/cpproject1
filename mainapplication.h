#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

#include <iostream>
#include "lib.h"
#include <ctime>
#include "vmcbf.h"
#include "vmcis.h"
#include <unistd.h>
#include "testdirectory.h"


class MainApplication
{
public:
    MainApplication(int _myrank, int _numprocs);
//    void calculateCFImportanceSampling();
//    void calculateClosedForm();
//    void alphabetavalues();
//    void steplength_manually();
//    void steplength_secant();
    void runApplication();

    int numprocs, myrank;

    void runSimulation();
    void simulateWithOutput(int nParticles, double alpha, double beta);
};

#endif // MAINAPPLICATION_H
