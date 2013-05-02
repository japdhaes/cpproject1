#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "wavefunction.h"
#include "slater.h"
#include "jastrow.h"

#include <armadillo>
using namespace std;
using namespace arma;

class Wavefunction
{
public:
    Wavefunction(int _nParticles, double _alpha, double _beta);

    void setAlpha(const double alpha);
    void setBeta(const double beta);
    void initialize(const mat &r);
    void acceptMove();
    void rejectMove();
    void setCurrentParticle(const int &i);

    double evaluate(const mat &r);
    double localEnergy(const mat &r);

    double getRatio();
    void setNewPos(const mat &r);
    double calcRatio();
protected:
    int nParticles;
    int nDimensions;
    int charge;
    double h, h2;
    int cp;

    mat rNew, rOld;
    Slater slater;
    Jastrow jastrow;
};

#endif // WAVEFUNCTION_H
