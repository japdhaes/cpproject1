#include "dimoleculewavefunction.h"

double alphafunction2(double alpha){
    double R = 0.5;
    double ans = 1+exp(-R/alpha)-alpha;
    return ans;
}

DimoleculeWavefunction::DimoleculeWavefunction(double _R, double _beta):
    alpha(0),
    beta(_beta),
    R(_R)
{
    alpha = rtsec(alphafunction2, 0, 100, 0.01);
}

double DimoleculeWavefunction::alphafunction(double alpha){
    double R = 0.5;
    double ans = 1+exp(-R/alpha)-alpha;
    return ans;
}
