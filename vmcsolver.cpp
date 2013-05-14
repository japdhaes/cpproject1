#include "vmcsolver.h"

//source github.com/ComputationalPhysics
VMCSolver::VMCSolver(int _myrank, int _numprocs, int _nParticles, double _alpha, double _beta):
    wf(_nParticles,_alpha,_beta),
    createOutput(false)
{
    this->nParticles=_nParticles;
    nDimensions=3;
    numprocs=_numprocs;
    myrank=_myrank;
    this->nCycles=1e5;
    idum=-1-myrank;
    h=1e-3;
    h2=1e6;
    local_nCycles=nCycles/numprocs;
}

VMCSolver::~VMCSolver(){

}


void VMCSolver::setAlphaBeta(double alpha, double beta){
    this->wf.setAlpha(alpha);
    this->wf.setBeta(beta);
}

double VMCSolver::runMonteCarloIntegration()
{
    nAccepted=0;
    nRejected=0;
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    waveFunctionOld = 0;
    waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

    // initialise trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j)=sqrt(h)*this->gaussianDeviate(&idum);
        }
    }

    initialize();
    // loop over Monte Carlo cycles
    int my_start=0;
    int my_end=local_nCycles;
    double r12 = 0;

    thermalize();

    for(int cycle = my_start; cycle < my_end; cycle++) {

//        if(cycle%(my_end/100)==0)
//            cout << "Currently in cycle "<< cycle<<endl;
        if(cycle%nrOfCyclesEachOutput==0 && createOutput && cycle >0){
            datalogger.flushData();
        }

        // Store the current value of the wave function
        waveFunctionOld = wf.evaluate(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            this->currentparticle=i;
            wf.setCurrentParticle(i);
            this->cycle(i);

        }
        //collect data
        deltaE = wf.localEnergyCF(rNew);
        if(createOutput){
            datalogger.addData(deltaE);
        }
        energySum += deltaE;
        energySquaredSum += deltaE*deltaE;
    }
    if(createOutput){
        datalogger.flushFinalData();
    }


    double energy = energySum/(local_nCycles);
    double totalenergy=0;
    int totalaccepted=0, totalrejected=0;

    //1 processor receives all information
    //MPI_Reduce(&energy, &totalenergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //all processors receive all information
    MPI_Allreduce(&energy, &totalenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nAccepted, &totalaccepted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nRejected, &totalrejected, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    double acceptratio=double(totalaccepted)/(totalrejected+totalaccepted);

    double energySquared = energySquaredSum/(local_nCycles);
    cout << "Energy: " << totalenergy/numprocs << " Energy (squared sum): " << energySquared << endl;
    cout << "Variance: "<< energySquared-pow(totalenergy/numprocs,2)<<endl;
    cout << "Total acceptratio: " <<acceptratio << endl;
    return totalenergy/numprocs;
}


double VMCSolver::gaussianDeviate(long *seed)
{
    double R, randomNormal;
    // Box-Muller transform
    R = sqrt(-2.0*log(ran2(seed)));
    randomNormal = R*cos(2.0*(pi) *ran2(seed));
    return randomNormal;
}

void VMCSolver::setCycles(const int &_nCycles){
    this->nCycles=_nCycles;
    local_nCycles=nCycles/numprocs;
}

void VMCSolver::setOutput(bool output){
    this->createOutput=output;
}

void VMCSolver::setOutputDirectory(const string &directoryName){
    this->outputDirectory=directoryName;
}

void VMCSolver::solverInitializer(){
    if(createOutput){
        ostringstream temp;
        temp << outputDirectory << "data"<< myrank<<".dat";
        this->outputFilename=temp.str();
        datalogger=Datalogger(temp.str(),nrOfCyclesEachOutput);
        datalogger.initialize();
    }
}

void VMCSolver::thermalize(){
    cout << "thermalizing "<<endl;
    for(int cycle = 0; cycle < thermalizingSteps; cycle++) {

//        if(cycle%(thermalizingSteps/50)==0)
//            cout << "Currently thermalizing in cycle "<< cycle<<endl;

        // Store the current value of the wave function
        waveFunctionOld = wf.evaluate(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            this->currentparticle=i;
            wf.setCurrentParticle(i);
            this->cycle(i);
        }
    }
}
