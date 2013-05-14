#ifndef Orbital_H
#define Orbital_H

#include <armadillo>

using namespace std;
using namespace arma;

class Orbital
{
public:
    Orbital(){};
    Orbital(double _alpha);
    void setAlpha(const double &newAlpha);
    double hydrogenWF(const vec &rvec, const int &qNum);
    vec3 gradient(const vec &rvec, const int &qNum);
    double laplacian(const vec &rvec, const int &qNum);
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
    double alphaGradient(const rowvec &rvec, const int &qNum);
private:
    int nDimensions;
    double wfCurrent;
    double alpha;
    double r, arg;
    vec3 grad;
    double lapl;
};

#endif // Orbital_H
