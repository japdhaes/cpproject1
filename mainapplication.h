#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

#include <iostream>
#include "lib.h"
#include <ctime>
#include "vmcsolver.h"
#include "hewf.h"
#include <unistd.h>
#include "bewf.h"

class MainApplication
{
public:
    MainApplication(int _myrank, int _numprocs);
    void calculateCFImportanceSampling();
    void calculateClosedForm();
    void alphabetavalues();
    void steplength_manually();
    void steplength_secant();
    void runApplication();

    int numprocs, myrank;


    void runBeryllium();
    void runBeryllium(double alpha, double beta);
};

#endif // MAINAPPLICATION_H
