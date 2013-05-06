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
}

double Wavefunction::evaluate(const mat &r)
{
    return /*jastrow.evaluate(r)**/slater.evaluate(r);
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
    double ans=this->slater.calcRatio()/**this->jastrow.calcRatio()*/;
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

double Wavefunction::localEnergyCF(const mat &r){
    double kinetic=localKineticCF();
    double potential=this->potentialEnergy(r);
//    cout << "CLOSED FORM: Kin " << kinetic << " pot " << potential << " tot " << kinetic+potential<<endl;
    return kinetic+potential;
}

double Wavefunction::localKineticCF(){
    //slide 173
    double kin=0.0;
//    cout << "SLATERDETERMINANT "<< slater.evaluate()<<endl;
    for(int i=0; i<this->nParticles;i++){
        kin+=slater.localLaplacian(i);
//        cout << "i "<<i<<" slaterlaplacian " <<slater.localLaplacian(i)<<endl;
//        kin+=jastrow.localLaplacian(i);
//        cout << "my laplacian "<<jastrow.localLaplacian(i)<<endl;
//        cout << "henrik laplacian "<<jastrow.getLaplaceRatio(rNew)<<endl;
////        cout << "i "<<i<<" slatergradient " << slater.localGradient(i)<<endl;
//        kin+=2.0*dot(slater.localGradient(i),jastrow.localGradient(i));
    }
    kin*=-0.5;
    return kin;
}

double Wavefunction::potentialEnergy(const mat &r){
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
//    double r12 = 0;
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = i + 1; j < nParticles; j++) {
//            r12 = 0;
//            for(int k = 0; k < nDimensions; k++) {
//                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
//            }
//            potentialEnergy += 1 / sqrt(r12);
//        }
//    }
    return potentialEnergy;
}

double Wavefunction::localEnergyNum(const mat &r)
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

mat Wavefunction::localGradient(){
    mat gradient=zeros(nParticles, nDimensions);
    gradient.zeros();
    for(int j=0; j<nParticles;j++){
        gradient.row(j) = slater.localGradient(j)+jastrow.localGradient(j);
    }
    return gradient;
}
