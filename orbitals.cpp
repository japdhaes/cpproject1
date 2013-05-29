#include "orbitals.h"

Orbitals::Orbitals():
    nDimensions(3),
    dphi(vec(nDimensions))
{
}

void Orbitals::setAlpha(const double &_alpha){
    this->alpha= _alpha;
}
