#ifndef DIMOLECULEWAVEFUNCTION_H
#define DIMOLECULEWAVEFUNCTION_H

#include "orbitals.h"
#include "lib.h"

class DimoleculeWavefunction
{
public:
    DimoleculeWavefunction(double _R, double _beta);
    double alphafunction(double alpha);
private:
    double alpha, beta, R;
};

#endif // DIMOLECULEWAVEFUNCTION_H
