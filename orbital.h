#ifndef Orbital_H
#define Orbital_H
#include "orbitals.h"
#include <armadillo>

using namespace std;
using namespace arma;

class Hydrogenic :public Orbitals
{
public:
    Hydrogenic():Orbitals(0){};
    Hydrogenic(double _alpha);
    virtual double wavefunction(const rowvec &rvec, const int &qNum) ;
    virtual void setAlpha(const double &_alpha);
    virtual rowvec gradient(const rowvec &rvec, const int &qNum);
    virtual double laplacian(const rowvec &rvec, const int &qNum) ;
    virtual double alphaGradient(const rowvec &rvec, const int &qNum);
    virtual void setR(const double &dist);
protected:
    double phi1s(const vec &rvec);
    double phi1s(const double &r);
    double phi2s(const vec &rvec);
    double phi2p(const vec &rvec, const int &k);
    vec3 dphi1s(const vec3 &rvec);
    vec3 dphi2s(const vec3 &rvec);
    vec3 dphi2p(const vec3 &rvec, const int &k);
    double ddphi1s(const vec &rvec);
    double ddphi2s(const vec &rvec);
    double ddphi2p(const vec &rvec, const int &k);
    double phi2s(const double &r);

};

#endif // Orbital_H
