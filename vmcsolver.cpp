#include "vmcsolver.h"
#include "zignor.c"
#include "zigrandom.c"

//source github.com/ComputationalPhysics
VMCSolver::VMCSolver(int _myrank, int _numprocs)
{
    stepLength=1.0;
    nDimensions=3;
    numprocs=_numprocs;
    myrank=_myrank;
    this->nCycles=1e6;
    idum=-1-myrank;
    h=1e-3;
    h2=1e6;
    local_nCycles=nCycles/numprocs;
    this->D=0.5;

    int temp = int(idum);
    RanNormalSetSeedZigVec(&temp, 100);
    nLocalTotalsteps=nParticles*local_nCycles;
}

VMCSolver::~VMCSolver(){

}

double VMCSolver::runMonteCarloIntegration()
{
    nAccepted=0;
    nRejected=0;
    //why a nparticles x ndimensions matrix? I have no clue what the point is of this structure.

    this->oldS.r = zeros<mat>(nParticles, nDimensions);
    this->oldS.rij = zeros<mat>(nParticles, nDimensions);
    this->oldS.wavefunction = 0.0;
    this->newS.r = zeros<mat>(nParticles, nDimensions);
    this->newS.rij = zeros<mat>(nParticles, nDimensions);
    this->newS.wavefunction = 0.0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            if(!importanceSampling){
                this->oldS.r(i,j) = stepLength * (ran2(&idum) - 0.5);
            } else {
//                rOld(i,j) = DRanNormalZigVec()*sqrt(h);
                this->oldS.r(i,j)=sqrt(h)*this->gaussianDeviate(&idum);
            }
        }
    }
    this->newS.r = this->oldS.r;

    if(importanceSampling){
        waveFunction(this->oldS);
//!!!        this->oldS.wavefunction = waveFunction(this->oldS.r);
//!!!        calculaterij(this->oldS.r, this->oldS.rij);
        calculaterij(this->oldS);
        this->newS.rij=this->oldS.rij;
//        fijOld = calculatefij(rijOld, )
        this->newS.qForce = this->oldS.qForce = quantumForce(this->oldS);
//!!!        qForceNew = qForceOld = quantumForce(rOld, waveFunctionOld) ;
    }



    // loop over Monte Carlo cycles
    int my_start=0;
    int my_end=local_nCycles;
    double r12 = 0;

    for(int cycle = my_start; cycle < my_end; cycle++) {
        // Store the current value of the wave function
        waveFunction(this->oldS);
//!!!        waveFunctionOld = waveFunction(rOld);


        // New position to test
        for(int i = 0; i < nParticles; i++) {

            if(importanceSampling){
                this->cycleIS(i);
            }
            else{
                this->cycleWithoutIS(i);
            }

            // update energies
            if(closedForm){
                deltaE = localEnergyClosedForm(this->newS);
//!!!                deltaE = localEnergyClosedForm(rNew);
            } else {
                deltaE = localEnergy(this->newS);
//!!!                deltaE = localEnergy(rNew);
            }

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;

        }
        //r12+=distij(rNew, 0, 1);
        r12+=norm(rNew.row(0),2);
    }

    r12/=local_nCycles;
    double total_r12=0;

    double energy = energySum/(local_nCycles * nParticles);
    double totalenergy=0;
    int totalaccepted=0, totalrejected=0;

    //1 processor receives all information
    //MPI_Reduce(&energy, &totalenergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //all processors receive all information
    MPI_Allreduce(&energy, &totalenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&r12, &total_r12, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nAccepted, &totalaccepted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nRejected, &totalrejected, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//    deltaP=double(totalaccepted)/(totalrejected+totalaccepted);

//    double energySquared = energySquaredSum/(local_nCycles * nParticles);
//    cout << "Energy: " << totalenergy/numprocs << " Energy (squared sum): " << energySquared << endl;
//    cout << "Optimal position electrons: " << total_r12/numprocs << endl;
//    cout << "Total deltaprobability: " <<deltaP << endl;
    return totalenergy/numprocs;
}

double VMCSolver::calcElementfij(const state &astate, int const &i, int const &j){
    mat rij=astate.rij;
    double a = ((i+j)%2==0)? 0.25 : 0.5;
    double b = rij(i,j);

    return a*b/(1.0+beta*b);
}


//double VMCSolver::calcElementfij(const mat &rij, int const i, int const j){
//    double a = ((i+j)%2==0)? 0.25 : 0.5;
//    double b = rij(i,j);

//    return a*b/(1.0+beta*b);
//}

void VMCSolver::calculatefij(state &astate){
    mat fij=astate.fij;
    mat rij=astate.rij;
    for(int i=0; i<nParticles;i++){
        for(int j=i+1; j<nParticles;j++){
            fij(i,j)=calcElementfij(astate,i,j);
            fij(j,i)=calcElementfij(astate,i,j);
        }
    }
    if(astate.fij(0,0)!=fij(0,0)){
        cout << "PROBLEM IN CALCFIJ"<<endl;
    }
}

//void VMCSolver::calculatefij(const mat &rij, mat &fij){
//    for(int i=0; i<nParticles;i++){
//        for(int j=i+1; j<nParticles;j++){
//            fij(i,j)=calcElementfij(rij,i,j);
//            fij(j,i)=calcElementfij(rij,i,j);
//        }
//    }
//}

void VMCSolver::cycleIS(const int &i){
    for(int j = 0; j < nDimensions; j++) {
//        rNew(i,j) = rOld(i,j) + DRanNormalZigVec()*sqrt(h)+qForceOld(i,j)*D*h;
        this->newS.r(i,j)=this->oldS.r(i,j)+ this->gaussianDeviate(&idum)*sqrt(h)+this->oldS.qForce(i,j)*D*h;
//        rNew(i,j) = rOld(i,j) + this->gaussianDeviate(&idum)*sqrt(h)+qForceOld(i,j)*D*h;
    }


    // Recalculate the value of the wave function
    this->newS.wavefunction=waveFunction(this->newS);
//    waveFunctionNew = waveFunction(rNew);
//!!!    qForceNew=this->quantumForce(rNew, waveFunctionNew);
    for(int j=0; j<nDimensions;j++){
        this->newS.qForce(i,j)=this->quantumForce(this->newS)(i,j);
    }

    // Check for step acceptance (if yes, update position, if no, reset position)
    //accepting
    double expratio = 0;
    for(int j=0; j<nDimensions;j++){
        expratio+=(this->oldS.qForce(i,j)+this->newS.qForce(i,j))*(D*h*(this->oldS.qForce(i,j)-this->newS.qForce(i,j))+2.0*(this->oldS.r(i,j)-this->newS.r(i,j)));
    }
    expratio/=4.0;
    expratio=exp(expratio);

    if(ran2(&idum) <= expratio*this->sdratio()*this->jastrowRatio(i)) {
        for(int j = 0; j < nDimensions; j++) {
            this->oldS.r(i,j) = this->newS.r(i,j);
            this->oldS.qForce(i,j)=this->newS.qForce(i,j);
            this->oldS.wavefunction = this->newS.wavefunction;
        }
        nAccepted++;
    }
    //rejecting
    else {
        for(int j = 0; j < nDimensions; j++) {
            this->newS.r(i,j) = this->oldS.r(i,j);
            this->newS.qForce(i,j)=this->oldS.qForce(i,j);
        }
        nRejected++;
    }
}

void VMCSolver::cycleWithoutIS(const int &i){
    for(int j = 0; j < nDimensions; j++) {
        this->newS.r(i,j) = this->oldS.r(i,j) + stepLength*(ran2(&idum) - 0.5);
    }

    // Recalculate the value of the wave function
    this->newS.wavefunction = waveFunction(this->newS);
//!!!    calculaterij(rNew, rijNew);
    calculaterij(this->newS);

    //accepting
    if(ran2(&idum) <= this->sdratio()*this->jastrowRatio(i)) {
        for(int j = 0; j < nDimensions; j++) {
            this->oldS.r(i,j) = this->newS.r(i,j);
            this->oldS.wavefunction = this->newS.wavefunction;
            this->oldS.rij(i,j)=this->newS.rij(i,j);
        }
        nAccepted++;
    }
    //rejecting
    else {
        for(int j = 0; j < nDimensions; j++) {
            this->newS.r(i,j) = this->oldS.r(i,j);
            this->newS.rij(i,j)=this->oldS.rij(i,j);
        }
        nRejected++;
    }
}

//void VMCSolver::cycleWithoutIS(const int &i){
//    for(int j = 0; j < nDimensions; j++) {
//        rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
//    }

//    // Recalculate the value of the wave function
//    waveFunctionNew = waveFunction(rNew);
//    calculaterij(rNew, rijNew);

//    //accepting
//    if(ran2(&idum) <= this->sdratio()*this->jastrowRatio(i)) {
//        for(int j = 0; j < nDimensions; j++) {
//            rOld(i,j) = rNew(i,j);
//            waveFunctionOld = waveFunctionNew;
//            rijOld(i,j)=rijNew(i,j);
//        }
//        nAccepted++;
//    }
//    //rejecting
//    else {
//        for(int j = 0; j < nDimensions; j++) {
//            rNew(i,j) = rOld(i,j);
//            rijNew(i,j)=rijOld(i,j);
//        }
//        nRejected++;
//    }
//}

mat VMCSolver::quantumForce(const state &astate){
    state rPlus;
    rPlus.r= zeros<mat>(nParticles, nDimensions);
    state rMinus;
    rMinus.r = zeros<mat>(nParticles, nDimensions);
    mat qforce(nParticles, nDimensions);
    double waveFunctionMinus, waveFunctionPlus;
    rPlus.r = rMinus.r =astate.r;

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus.r(i,j) += h;
            rMinus.r(i,j) -= h;
            waveFunctionMinus = waveFunction(rMinus);
            waveFunctionPlus = waveFunction(rPlus);
            qforce(i,j) = waveFunctionPlus-waveFunctionMinus;
            rPlus.r(i,j) = astate.r(i,j);
            rMinus.r(i,j) = astate.r(i,j);
        }
    }
    qforce/=(astate.wavefunction*h);
    return qforce;
}


//mat VMCSolver::quantumForce(const mat &r, const double &wavefunction){
//    mat qforce(nParticles, nDimensions);
//    mat rPlus = zeros<mat>(nParticles, nDimensions);
//    mat rMinus = zeros<mat>(nParticles, nDimensions);
//    double waveFunctionMinus, waveFunctionPlus;
//    rPlus = rMinus =r;

//    for(int i = 0; i < nParticles; i++) {
//        for(int j = 0; j < nDimensions; j++) {
//            rPlus(i,j) += h;
//            rMinus(i,j) -= h;
//            waveFunctionMinus = waveFunction(rMinus);
//            waveFunctionPlus = waveFunction(rPlus);
//            qforce(i,j) = waveFunctionPlus-waveFunctionMinus;
//            rPlus(i,j) = r(i,j);
//            rMinus(i,j) = r(i,j);
//        }
//    }
//    qforce/=(wavefunction*h);
//    return qforce;
//}

double VMCSolver::distij(const state &astate, const int &i, const int &j){
    const mat r = astate.r;
    return norm(r.row(i)-r.row(j),2);
}

//double VMCSolver::distij(const mat &r, const int i, const int j){
//    return norm(r.row(i)-r.row(j),2);
//}

void VMCSolver::updaterij(state &astate, const int &j){
    mat rij = astate.rij;
    const mat r = astate.r;
    for(int i=0; i<j;i++){
        rij(i,j)=0.0;
        for(int k=0; k<nDimensions;k++){
            rij(i,j)+=(r(i,k)-r(j,k))*(r(i,k)-r(j,k));
        }
        rij(i,j)=sqrt(rij(i,j));

    }
    for(int i=j+1; i<nParticles;i++){
        rij(j,i)=0.0;
        for(int k=0; k<nDimensions;k++){
            rij(j,i)+=(r(i,k)-r(j,k))*(r(i,k)-r(j,k));
        }
        rij(j,i)=sqrt(rij(j,i));
    }
}

//void VMCSolver::updaterij(const mat &r, mat &rij, const int j){

//    for(int i=0; i<j;i++){
//        rij(i,j)=0.0;
//        for(int k=0; k<nDimensions;k++){
//            rij(i,j)+=(r(i,k)-r(j,k))*(r(i,k)-r(j,k));
//        }
//        rij(i,j)=sqrt(rij(i,j));

//    }
//    for(int i=j+1; i<nParticles;i++){
//        rij(j,i)=0.0;
//        for(int k=0; k<nDimensions;k++){
//            rij(j,i)+=(r(i,k)-r(j,k))*(r(i,k)-r(j,k));
//        }
//        rij(j,i)=sqrt(rij(j,i));
//    }
//}

void VMCSolver::calculaterij(state &astate){
    mat rij = astate.rij;
    const mat r = astate.r;
    rij.zeros();
    for(int i=0; i<nParticles; i++){
        for(int j=i+1; j<nParticles;j++){
            rij(i,j)=0.0;
            for(int k=0; k<nDimensions;k++){
                rij(i,j)+=(r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            rij(i,j)=sqrt(rij(i,j));
        }
    }
}

//void VMCSolver::calculaterij(const mat &r, mat &rij){
//    rij.zeros();
//    for(int i=0; i<nParticles; i++){
//        for(int j=i+1; j<nParticles;j++){
//            rij(i,j)=0.0;
//            for(int k=0; k<nDimensions;k++){
//                rij(i,j)+=(r(i,k)-r(j,k))*(r(i,k)-r(j,k));
//            }
//            rij(i,j)=sqrt(rij(i,j));
//        }
//    }
//}


//void VMCSolver::calculaterij(const mat &r, mat &rij){
//    rij.zeros();
//    for(int i=0; i<nParticles; i++){
//        for(int j=i+1; j<nParticles;j++){
//            rij(i,j)=0.0;
//            for(int k=0; k<nDimensions;k++){
//                rij(i,j)+=(r(i,k)-r(j,k))*(r(i,k)-r(j,k));
//            }
//            rij(i,j)=sqrt(rij(i,j));
//        }
//    }
//}


double VMCSolver::localEnergy(const state &astate)
{
    state rPlus, rMinus;
    rPlus.r = zeros<mat>(nParticles, nDimensions);
    rMinus.r = zeros<mat>(nParticles, nDimensions);
    const mat r = astate.r;

    rPlus.r = rMinus.r = astate.r;

    rMinus.wavefunction = 0;
    rPlus.wavefunction = 0;

    double waveFunctionCurrent = waveFunction(astate);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus.r(i,j) += h;
            rMinus.r(i,j) -= h;
            rMinus.wavefunction = waveFunction(rMinus);
            rPlus.wavefunction = waveFunction(rPlus);
            kineticEnergy -= (rMinus.wavefunction + rPlus.wavefunction - 2 * waveFunctionCurrent);
            rPlus.r(i,j) = astate.r(i,j);
            rMinus.r(i,j) = astate.r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += astate.r(i,j)*astate.r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (astate.r(i,k) - astate.r(j,k)) * (astate.r(i,k) - astate.r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    //cout << "kinetic energy "<<kineticEnergy<< " potential energy " <<potentialEnergy<<endl;
    return kineticEnergy + potentialEnergy;
}

double VMCSolver::gaussianDeviate(long *seed)
{
    double R, randomNormal;
    // Box-Muller transform
//    randomUniform << ran2(seed) << ran2(seed);
//    R = sqrt(-2*log(randomUniform(0)));
//    randomNormal(0) = R*cos(2*pi*randomUniform(1));
//    randomNormal(1) = R*sin(2*pi*randomUniform(1))

    R = sqrt(-2.0*log(ran2(seed)));
    randomNormal = R*cos(2.0*(pi) *ran2(seed));
    return randomNormal;
}

double VMCSolver::hydrogenWF(const int &i, vec3 &r)
{
    double wf;
    switch (i)
    {
    case 0:
        wf = phi1s(r);
        break;
    case 1:
        wf = phi2s(r);
        break;
    default:
        cout << "wrong implementation of hydrogenWF"<<endl;
    }
    return wf;

}

double VMCSolver::phi1s(const double position){
    return exp(-alpha*position);
}

double VMCSolver::phi2s(const double position){
    double arg=alpha*position*0.5;
    return (1.0-arg)*exp(-arg);
}
double VMCSolver::phi1s(const vec3 &position){
    double r1=0;
    for(int i=0; i<nDimensions;i++){
        r1+=(position(i))*position(i);
    }
    r1=sqrt(r1);
    return exp(-alpha*r1);
}

double VMCSolver::phi2s(const vec3 &position){
    double arg=0;
    for(int i=0; i<nDimensions;i++){
        arg+=(position(i))*position(i);
    }
    arg=sqrt(arg);
    arg=alpha*arg*0.5;
    return (1-arg)*exp(-arg);
}
