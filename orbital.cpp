#include "orbital.h"

Hydrogenic::Hydrogenic(double _alpha)
{
    this->alpha=_alpha;
    this->nDimensions=3;
}

void Hydrogenic::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
}

double Hydrogenic::wavefunction(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0 :
        wfCurrent = phi1s(rvec);
        break;
    case 1 :
        wfCurrent = phi2s(rvec);
        break;
    case 2:
        wfCurrent = phi2p(rvec, 0);
        break;
    case 3:
        wfCurrent = phi2p(rvec, 1);
        break;
    case 4:
        wfCurrent = phi2p(rvec, 2);
        break;
    default :
        // Process for all other cases.
        cout << "! We don't have this orbital yet!" << endl;
        exit(1);
    }

    return wfCurrent;
}

rowvec Hydrogenic::gradient(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0:
        grad = dphi1s(rvec);
        break;
    case 1:
        grad = dphi2s(rvec);
        break;
    case 2:
        grad = dphi2p(rvec, 0);
        break;
    case 3:
        grad = dphi2p(rvec, 1);
        break;
    case 4:
        grad = dphi2p(rvec, 2);
        break;
    default:
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }

    return grad;
}

double Hydrogenic::laplacian(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0:
        lapl = ddphi1s(rvec);
        break;
    case 1:
        lapl = ddphi2s(rvec);
        break;
    case 2:
        lapl = ddphi2p(rvec, 0);
        break;
    case 3:
        lapl = ddphi2p(rvec, 1);
        break;
    case 4:
        lapl = ddphi2p(rvec, 2);
        break;
    default:
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }

    return lapl;
}

double Hydrogenic::phi1s(const vec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    return exp(-alpha*r);
}

double Hydrogenic::phi1s(const double &rr)
{
    return exp(-alpha*rr);
}

double Hydrogenic::phi2s(const vec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);
    arg = alpha*r*0.5;

    return (1.0 - arg)*exp(-arg);
}

double Hydrogenic::phi2s(const double &r)
{
    double arg = alpha*r*0.5;
    return (1.0 - arg)*exp(-arg);
}

double Hydrogenic::phi2p(const vec &rvec, const int &k)
{
    r = 0.0;
    for (int i = 0; i<nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return rvec(k)*exp(-0.5*alpha*r);
}

vec3 Hydrogenic::dphi1s(const vec3 &rvec)
{
    vec3 dphi;
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    dphi = (-alpha/r*exp(-alpha*r))*rvec;

    return dphi;
}

vec3 Hydrogenic::dphi2s(const vec3 &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return alpha*rvec*(alpha*r - 4.0)*exp(-alpha*r/2.0)/(4.0*r);
}

vec3 Hydrogenic::dphi2p(const vec3 &rvec, const int &k)
{
    vec3 dphi;
    double r2 = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r2 += rvec(i)*rvec(i);
    }
    r = sqrt(r2);
    dphi = -alpha * rvec(k) * rvec;
//    for (int i = 0; i < nDimensions; i++)
//    {
//        dphi(i) = -alpha*rvec(i)*rvec(k);
//    }
    dphi(k) += 2*r;
    dphi *= exp(-alpha*r/2.0)/(2*r);
    return dphi;
}

double Hydrogenic::ddphi1s(const vec &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    return alpha*(alpha*r - 2.0)*exp(-alpha*r)/r;
}

double Hydrogenic::ddphi2s(const vec &rvec)
{
    double r2 = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r2 += rvec(i)*rvec(i);
    }
    r = sqrt(r2);

    return -alpha*(alpha*alpha*r2 - 10*alpha*r + 16)*exp(-alpha*r/2)/(8*r);
}

double Hydrogenic::ddphi2p(const vec &rvec, const int &k)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return alpha*rvec(k)*(alpha*r - 8.0)*exp(-alpha*r/2.0)/(4.0*r);
}

double Hydrogenic::alphaGradient(const rowvec &rvec, const int &qNum)
{
    double r = 0.0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    switch (qNum)
    {
    case 0 :
        return -r*exp(-alpha*r);
    case 1 :
        return r*(alpha*r - 4.0)*exp(-alpha*r/2.0)/4.0;
    case 2 :
        return -rvec(0)*(alpha*r - 2.0)*exp(-alpha*r/2.0)/2.0;
    case 3 :
        return -rvec(1)*(alpha*r - 2.0)*exp(-alpha*r/2.0)/2.0;
    case 4 :
        return -rvec(2)*(alpha*r - 2.0)*exp(-alpha*r/2.0)/2.0;
    default :
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }
}

void Hydrogenic::setR(const double &dist)
{
}

