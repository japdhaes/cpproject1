#ifndef DIMOLECULEWAVEFUNCTION_H
#define DIMOLECULEWAVEFUNCTION_H

#include "orbitals.h"
#include "lib.h"
#include "slater.h"
#include "jastrow.h"

class DimoleculeWavefunction
{
public:
    DimoleculeWavefunction(int _nParticles, double alpha, double _beta);
    double alphafunction(double alpha);
    void setAlpha(const double &alpha);
    void setBeta(const double beta);
    void initialize(const mat &r);
    void acceptMove();
    void rejectMove();
    void setCurrentParticle(const int &i);

    double evaluate(const mat &r);

    double localEnergyNum(const mat &r);

    double getRatio();
    void setNewPos(const mat &r);
    double calcRatio();

    double localEnergyCF(const mat &r);
    double localKineticCF();
    double potentialEnergy(const mat &r);
    mat localGradient();
    void setR(double _R);
protected:
    double alpha, beta, R;
    int nParticles;
    int nDimensions;
    int charge;
    double h, h2;
    int cp;

    Slater slater;
    Jastrow jastrow;
};

#endif // DIMOLECULEWAVEFUNCTION_H
