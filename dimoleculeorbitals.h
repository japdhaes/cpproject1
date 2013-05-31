#ifndef DIMOLECULEORBITALS_H
#define DIMOLECULEORBITALS_H

#include <armadillo>
#include <orbital.h>
using namespace std;
using namespace arma;

class DimoleculeOrbitals : public Orbitals
{
public:
    DimoleculeOrbitals(double _alpha, double _R);
    void setAlpha(const double &newAlpha);
    void setR(const double &dist);
    double wavefunction(const rowvec &rvec, const int &qNum);
    rowvec gradient(const rowvec &rvec, const int &qNum);
    double laplacian(const rowvec &rvec, const int &qNum);
    double alphaGradient(const rowvec &rvec, const int &qNum);

private:
    Hydrogenic hydrogenic;
    rowvec R;
};

#endif // DIMOLECULEORBITALS_H
