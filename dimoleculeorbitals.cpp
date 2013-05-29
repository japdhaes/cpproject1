#include "dimoleculeorbitals.h"

DimoleculeOrbitals::DimoleculeOrbitals(double _alpha, double _dist):
    Orbitals(),
    R(zeros<rowvec>(nDimensions)),
    alpha(_alpha),
    hydrogenic(Hydrogenic(_alpha))
{
    R(0)=_dist/2;
}

void DimoleculeOrbitals::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
    hydrogenic.setAlpha(newAlpha);
}

void DimoleculeOrbitals::setR(const double &dist)
{
    R(0) = dist/2.0;
}

double DimoleculeOrbitals::wavefunction(const rowvec &rvec, const int &qNum)
{
    return hydrogenic.wavefunction(rvec + R, qNum)
            + hydrogenic.wavefunction(rvec - R, qNum);
}

rowvec DimoleculeOrbitals::gradient(const rowvec &rvec, const int &qNum)
{
    return hydrogenic.gradient(rvec + R, qNum)
            + hydrogenic.gradient(rvec - R, qNum);
}

double DimoleculeOrbitals::laplacian(const rowvec &rvec, const int &qNum)
{
    return hydrogenic.laplacian(rvec + R, qNum)
            + hydrogenic.laplacian(rvec - R, qNum);
}

double DimoleculeOrbitals::alphaGradient(const rowvec &rvec, const int &qNum)
{
    return hydrogenic.alphaGradient(rvec + R, qNum)
            + hydrogenic.alphaGradient(rvec - R, qNum);
}
