#include "wavefunction.h"

Wavefunction::Wavefunction(int _nParticles, double _alpha, double _beta):
    slater(_alpha, _nParticles),
    jastrow(_beta, _nParticles),
    charge(_nParticles),
    nParticles(_nParticles),
    nDimensions(3),
    h(1e-3),
    h2(1e6)
{
}

void Wavefunction::setAlpha(const double alpha){
    this->slater.setAlpha(alpha);
}

void Wavefunction::setBeta(const double beta){
    this->jastrow.setBeta(beta);
}

void Wavefunction::initialize(const mat &r){
    this->slater.initialize(r);
    this->jastrow.initialize(r);
    cout << sizeof(slater) << " "<<sizeof(jastrow)<<endl;
}

double Wavefunction::evaluate(const mat &r)
{
    return jastrow.evaluate(r)*slater.evaluate(r);
}

double Wavefunction::localEnergy(const mat &r)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = evaluate(r);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = evaluate(rMinus);
            waveFunctionPlus = evaluate(rPlus);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    cout << "kinetic energy "<<kineticEnergy<< " potential energy " <<potentialEnergy<<" total " << kineticEnergy + potentialEnergy<<endl;
    return kineticEnergy + potentialEnergy;
}


void Wavefunction::setCurrentParticle(const int &i)
{
    this->cp=i;
    this->slater.setCurrentParticle(i);
    this->jastrow.setCurrentParticle(i);
}

void Wavefunction::setNewPos(const mat &r){
    this->slater.setNewPosition(r);
    this->jastrow.setNewPosition(r);
}

double Wavefunction::calcRatio(){
    double ans=this->slater.calcRatio()*this->jastrow.calcRatio();
    return ans*ans;
}

double Wavefunction::getRatio(){
    double ans=this->slater.getRatio()*this->jastrow.getRatio();
    return ans*ans;
}

void Wavefunction::acceptMove(){
    this->slater.acceptMove();
    this->jastrow.acceptMove();
}

void Wavefunction::rejectMove(){
    this->slater.rejectMove();
    this->jastrow.rejectMove();
}
