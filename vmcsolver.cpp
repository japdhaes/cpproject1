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

    nLocalTotalsteps=nParticles*local_nCycles;
}

VMCSolver::~VMCSolver(){

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

    for(int cycle = my_start; cycle < my_end; cycle++) {

        if(cycle%(my_end/100)==0)
            cout << "Currently in cycle "<< cycle<<endl;
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

            //collect data
            deltaE = wf.localEnergyCF(rNew);
            if(createOutput){
                datalogger.addData(deltaE);
            }
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }

        //r12+=distij(rNew, 0, 1);
//        r12+=norm(rNew.row(0),2);
    }
    if(createOutput){
        datalogger.flushFinalData();
    }

//    r12/=local_nCycles;
//    double total_r12=0;

    double energy = energySum/(local_nCycles * nParticles);
    double totalenergy=0;
    int totalaccepted=0, totalrejected=0;

    //1 processor receives all information
    //MPI_Reduce(&energy, &totalenergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //all processors receive all information
    MPI_Allreduce(&energy, &totalenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&r12, &total_r12, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nAccepted, &totalaccepted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nRejected, &totalrejected, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    double acceptratio=double(totalaccepted)/(totalrejected+totalaccepted);

    double energySquared = energySquaredSum/(local_nCycles * nParticles);
    cout << "Energy: " << totalenergy/numprocs << " Energy (squared sum): " << energySquared << endl;
    cout << "Variance: "<< energySquared-pow(totalenergy/numprocs,2)<<endl;
//    cout << "Optimal position electrons: " << total_r12/numprocs << endl;
    cout << "Total acceptratio: " <<acceptratio << endl;
//    double totalenergy=0.0;
    return totalenergy/numprocs;
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

void VMCSolver::setCycles(const int &_nCycles){
    this->nCycles=_nCycles;
    local_nCycles=nCycles/numprocs;
    cout << local_nCycles<<endl;
    nLocalTotalsteps=nParticles*local_nCycles;
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
        datalogger=Datalogger(temp.str(),nrOfCyclesEachOutput*nParticles);
        datalogger.initialize();
    }
}
