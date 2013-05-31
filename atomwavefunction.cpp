#include "atomwavefunction.h"

AtomWavefunction::AtomWavefunction(int _nParticles, double _alpha, double _beta):
    slater(_alpha, _nParticles, 0),
    jastrow(_beta, _nParticles),
    charge(_nParticles),
    nParticles(_nParticles),
    nDimensions(3),
    h(5e-3),
    h2(1.0/h/h)
{
}

void AtomWavefunction::setAlpha(const double alpha){
    this->slater.setAlpha(alpha);
}

void AtomWavefunction::setR(double R){}

void AtomWavefunction::setBeta(const double beta){
    if(beta<0){
        cout << "Setting a negative beta value. Exiting."<<endl;
        exit(1);
    }
    this->jastrow.setBeta(beta);
}

void AtomWavefunction::initialize(const mat &r){
    this->slater.initialize(r);
    this->jastrow.initialize(r);
}

double AtomWavefunction::evaluate(const mat &r)
{
    return jastrow.evaluate(r)*slater.evaluate(r);
}



void AtomWavefunction::setCurrentParticle(const int &i)
{
    this->cp=i;
    this->slater.setCurrentParticle(i);
    this->jastrow.setCurrentParticle(i);
}

void AtomWavefunction::setNewPos(const mat &r){
    this->slater.setNewPosition(r);
    this->jastrow.setNewPosition(r);
}

double AtomWavefunction::calcRatio(){
    double ans=this->slater.calcRatio()*this->jastrow.calcRatio();
    return ans*ans;
}

double AtomWavefunction::getRatio(){
    double ans=this->slater.getRatio()*this->jastrow.getRatio();
    return ans*ans;
}

void AtomWavefunction::acceptMove(){
    this->slater.acceptMove();
    this->jastrow.acceptMove();
}

void AtomWavefunction::rejectMove(){
    this->slater.rejectMove();
    this->jastrow.rejectMove();
}

double AtomWavefunction::localEnergyCF(const mat &r){
    double kinetic=localKineticCF();
    double potential=this->potentialEnergy(r);
//    cout << r << endl;
//    cout << "CLOSED FORM: Kin " << kinetic << " pot " << potential << " tot " << kinetic+potential<<endl;
    return kinetic+potential;
}

double AtomWavefunction::localKineticCF(){
    //slide 173
    double kin=0.0;
    for(int i=0; i<this->nParticles;i++){
        kin+=slater.localLaplacian(i);
        kin+=jastrow.localLaplacian(i);
        kin+=2.0*dot(slater.localGradient(i),jastrow.localGradient(i));
    }
    kin*=-0.5;
    return kin;
}

double AtomWavefunction::potentialEnergy(const mat &r){
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    // Contribution from electron-nucleus potential
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
    return potentialEnergy;
}

double AtomWavefunction::localEnergyNum(const mat &r)
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
    double potentialEnergy = this->potentialEnergy(r);

    //cout << "kinetic energy "<<kineticEnergy<< " potential energy " <<potentialEnergy<<" total " << kineticEnergy + potentialEnergy<<endl;
    return kineticEnergy + potentialEnergy;
}

mat AtomWavefunction::localGradient(){
    mat gradient=zeros(nParticles, nDimensions);
    gradient.zeros();
    for(int j=0; j<nParticles;j++){
        gradient.row(j) = slater.localGradient(j)+jastrow.localGradient(j);
    }
    return gradient;
}
