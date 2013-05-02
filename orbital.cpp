#include "orbital.h"

Orbital::Orbital(double _alpha)
{
    this->alpha=_alpha;
    this->nDimensions=3;
}

void Orbital::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
}

double Orbital::hydrogenWF(const vec &rvec, const int &qNum)
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

vec3 Orbital::gradient(const vec &rvec, const int &qNum)
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

double Orbital::laplacian(const vec &rvec, const int &qNum)
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

double Orbital::phi1s(const vec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    return exp(-alpha*r);
}

double Orbital::phi1s(const double &rr)
{
    return exp(-alpha*rr);
}

double Orbital::phi2s(const vec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);
    arg = alpha*r*0.5;

    return (1.0 - arg)*exp(-arg);
}

double Orbital::phi2s(const double &r)
{
    double arg = alpha*r*0.5;
    return (1.0 - arg)*exp(-arg);
}

double Orbital::phi2p(const vec &rvec, const int &k)
{
    r = 0.0;
    for (int i = 0; i<nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return alpha*rvec(k)*exp(-0.5*alpha*r);
}

vec3 Orbital::dphi1s(const vec3 &rvec)
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

vec3 Orbital::dphi2s(const vec3 &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return alpha*rvec*(alpha*r - 4.0)*exp(-alpha*r/2.0)/(4.0*r);
}

vec3 Orbital::dphi2p(const vec3 &rvec, const int &k)
{
    vec3 dphi;
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    double arg = exp(-alpha*r/2.0);
    for (int i = 0; i < nDimensions; i++)
    {
        dphi = -alpha*alpha/(2.0*r)*arg*rvec(i)*rvec(k);
    }
    dphi(k) += alpha*arg;

    return dphi;
}

double Orbital::ddphi1s(const vec &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    return alpha*(alpha*r - 2.0)*exp(-alpha*r)/r;
}

double Orbital::ddphi2s(const vec &rvec)
{
    double r2 = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r2 += rvec(i)*rvec(i);
    }
    r = sqrt(r2);

    return -alpha*(alpha*alpha*r2 - 10*alpha*r + 16)*exp(-alpha*r/2)/(8*r);
}

double Orbital::ddphi2p(const vec &rvec, const int &k)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return pow(alpha, 2.0)*rvec(k)*(alpha*r - 8.0)*exp(-alpha*r/2.0)/(4.0*r);
}

