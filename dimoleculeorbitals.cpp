#include "dimoleculeorbitals.h"

DimoleculeOrbitals::DimoleculeOrbitals(double _alpha, double _dist):
    Orbitals(_alpha),
    R(zeros<rowvec>(nDimensions)),
    hydrogenic(Hydrogenic(_alpha))
{
    R(0)=_dist/2;
}

void DimoleculeOrbitals::setR(const double &_R)
{
    R(0) = _R;
}
void DimoleculeOrbitals::setAlpha(const double &_alpha){
    this->alpha= _alpha;
    hydrogenic.setAlpha(_alpha);
}

double DimoleculeOrbitals::wavefunction(const rowvec &rvec, const int &qNum)
{
    if (qNum%2 == 0){
        return hydrogenic.wavefunction(rvec + R, qNum/2)
                + hydrogenic.wavefunction(rvec - R, qNum/2);
    }
    else{
        return hydrogenic.wavefunction(rvec + R, qNum/2)
                - hydrogenic.wavefunction(rvec - R, qNum/2);
    }
}

rowvec DimoleculeOrbitals::gradient(const rowvec &rvec, const int &qNum)
{
    if (qNum%2 == 0){
        return hydrogenic.gradient(rvec + R, qNum/2)
            + hydrogenic.gradient(rvec - R, qNum/2);
    }
    else{
        return hydrogenic.gradient(rvec + R, qNum/2)
            - hydrogenic.gradient(rvec - R, qNum/2);
    }
}

double DimoleculeOrbitals::laplacian(const rowvec &rvec, const int &qNum)
{
    if (qNum%2 == 0){
        return hydrogenic.laplacian(rvec + R, qNum/2)
            + hydrogenic.laplacian(rvec - R, qNum/2);
    }
    else{
        return hydrogenic.laplacian(rvec + R, qNum/2)
            - hydrogenic.laplacian(rvec - R, qNum/2);
    }
}

double DimoleculeOrbitals::alphaGradient(const rowvec &rvec, const int &qNum)
{
    return hydrogenic.alphaGradient(rvec + R, qNum)
            + hydrogenic.alphaGradient(rvec - R, qNum);
}
