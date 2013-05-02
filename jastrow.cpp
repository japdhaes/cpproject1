#include "jastrow.h"

Jastrow::Jastrow()
{
}

Jastrow::Jastrow(double _beta, int _nParticles)
{
    this->beta=_beta;
    this->nParticles=_nParticles;
    this->nDimensions=3;
}

void Jastrow::setBeta(double _beta)
{
    this->beta=_beta;
}

void Jastrow::initialize(const mat &r)
{
}

double Jastrow::evaluate(const mat &r)
{
    return 1.0;
}

void Jastrow::setCurrentParticle(const int &i)
{
    this->cp=i;
}

void Jastrow::setNewPosition(const mat &r){

}

double Jastrow::calcRatio()
{
    return 1.0;
}


double Jastrow::getRatio(){
    return 1.0;
}

void Jastrow::acceptMove()
{
}

void Jastrow::rejectMove()
{
}
