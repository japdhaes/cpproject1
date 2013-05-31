#include "dimoleculewavefunction.h"

double alphafunction2(double alpha, double R){
    double ans = 1+exp(-R/alpha)-alpha;
    return ans;
}

DimoleculeWavefunction::DimoleculeWavefunction(int _nParticles, double _alpha, double _beta):
    alpha(_alpha),
    beta(_beta),
    R(0.0),
    slater(_alpha, _nParticles, 1),
    jastrow(_beta, _nParticles),
    charge(_nParticles/2),
    nParticles(_nParticles),
    nDimensions(3),
    h(5e-3),
    h2(1.0/h/h)

{

}

void DimoleculeWavefunction::setR(double _R){
    this->R=_R;
    slater.setR(_R);
}

double DimoleculeWavefunction::alphafunction(double alpha){
    double R = 0.5;
    double ans = 1+exp(-R/alpha)-alpha;
    return ans;
}

void DimoleculeWavefunction::setAlpha(const double &alpha){
    this->alpha=alpha;
    this->slater.setAlpha(alpha);
}

void DimoleculeWavefunction::setBeta(const double beta){
    if(beta<0){
        cout << "Setting a negative beta value. Exiting."<<endl;
        exit(1);
    }
    this->jastrow.setBeta(beta);
}

void DimoleculeWavefunction::initialize(const mat &r){
    this->slater.initialize(r);
    this->jastrow.initialize(r);
}

double DimoleculeWavefunction::evaluate(const mat &r)
{
    return jastrow.evaluate(r)*slater.evaluate(r);
}



void DimoleculeWavefunction::setCurrentParticle(const int &i)
{
    this->cp=i;
    this->slater.setCurrentParticle(i);
    this->jastrow.setCurrentParticle(i);
}

void DimoleculeWavefunction::setNewPos(const mat &r){
    this->slater.setNewPosition(r);
    this->jastrow.setNewPosition(r);
}

double DimoleculeWavefunction::calcRatio(){
    double ans=this->slater.calcRatio()*this->jastrow.calcRatio();
    return ans*ans;
}

double DimoleculeWavefunction::getRatio(){
    double ans=this->slater.getRatio()*this->jastrow.getRatio();
    return ans*ans;
}

void DimoleculeWavefunction::acceptMove(){
    this->slater.acceptMove();
    this->jastrow.acceptMove();
}

void DimoleculeWavefunction::rejectMove(){
    this->slater.rejectMove();
    this->jastrow.rejectMove();
}

double DimoleculeWavefunction::localEnergyCF(const mat &r){
    double kinetic=localKineticCF();
    double potential=this->potentialEnergy(r);
//    cout << "CLOSED FORM: Kin " << kinetic << " pot " << potential << " tot " << kinetic+potential<<endl;
//    cout << evaluate(r )<<endl;
    return kinetic+potential;
}

double DimoleculeWavefunction::localKineticCF(){
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

double DimoleculeWavefunction::potentialEnergy(const mat &r){
    //cout << "N="<<nParticles<<" Z="<<charge<<" alpha="<<alpha<<" beta="<<beta<<" R="<<R<<endl;
    //cout << r<<endl;
    double potentialEnergy = 0;
    vec3 rNucleus; rNucleus << R << 0 << 0 ;
    // Contribution from electron-nucleus potential
    for(int i = 0; i < nParticles; i++) {
        double rp1 = 0;
        double rp2 = 0;
        for(int j = 0; j < nDimensions; j++) {
            rp1 += (r(i,j) + rNucleus(j))*(r(i,j) + rNucleus(j));
            rp2 += (r(i,j) - rNucleus(j))*(r(i,j) - rNucleus(j));
        }
        potentialEnergy -= charge*(1.0/sqrt(rp1) + 1.0/sqrt(rp2));
    }
    double currentpot = potentialEnergy;
    //cout << "eN potential=" <<currentpot<<endl;

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
    //cout << "ee potential=" <<potentialEnergy-currentpot<<endl;
    // Contribution from nucleus-nucleus potential
    //cout << "NN potential="<<charge*charge/(2.0*R)<<endl;
    potentialEnergy += charge*charge/(2.0*R);
    return potentialEnergy;
}

double DimoleculeWavefunction::localEnergyNum(const mat &r)
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

mat DimoleculeWavefunction::localGradient(){
    mat gradient=zeros(nParticles, nDimensions);
    gradient.zeros();
    for(int j=0; j<nParticles;j++){
        gradient.row(j) = slater.localGradient(j)+jastrow.localGradient(j);
    }
    return gradient;
}
