#include "vmcsolver.h"

//source github.com/ComputationalPhysics
VMCSolver::VMCSolver(int _myrank, int _numprocs, int _nParticles, double _alpha, double _beta):
    wf(_nParticles,_alpha,_beta)
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
        if(cycle%100==0)
            cout << "Currently in cycle "<< cycle<<endl;
        // Store the current value of the wave function
        waveFunctionOld = wf.evaluate(rOld);


        // New position to test
        for(int i = 0; i < nParticles; i++) {
            this->currentparticle=i;
            wf.setCurrentParticle(i);
            this->cycle(i);

            // update energies
//            if(closedForm){
//                deltaE = localEnergyClosedForm(rNew);
//            } else {
//            cout << rNew<<endl;
                deltaE = wf.localEnergy(rNew);
//                cout << deltaE<<endl;
//            }

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;

        }
        //r12+=distij(rNew, 0, 1);
//        r12+=norm(rNew.row(0),2);
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
//    deltaP=double(totalaccepted)/(totalrejected+totalaccepted);

    double energySquared = energySquaredSum/(local_nCycles * nParticles);
    cout << "Energy: " << totalenergy/numprocs << " Energy (squared sum): " << energySquared << endl;
//    cout << "Optimal position electrons: " << total_r12/numprocs << endl;
//    cout << "Total deltaprobability: " <<deltaP << endl;
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
