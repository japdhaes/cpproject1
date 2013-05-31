#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>

using namespace std;
using namespace arma;

class Orbitals
{
public:
    Orbitals(double _alpha);
    virtual void setAlpha(const double &_alpha)=0;
    virtual void setR(const double &dist) = 0;
    virtual double wavefunction(const rowvec &rvec, const int &qNum)=0 ;
    virtual rowvec gradient(const rowvec &rvec, const int &qNum) =0;
    virtual double laplacian(const rowvec &rvec, const int &qNum) =0;
    virtual double alphaGradient(const rowvec &rvec, const int &qNum) =0;
protected:
    int nDimensions;
    double wfCurrent;
    double alpha;
    double r, arg;
    rowvec grad;
    double lapl;
    vec dphi;

};

#endif // ORBITALS_H

#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>

using namespace std;
using namespace arma;

class Orbitals
{
public:



};

#endif // ORBITALS_H
